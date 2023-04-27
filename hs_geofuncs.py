import numpy as np
import pandas as pd
import time
import pystac_client
import planetary_computer
import rioxarray
import matplotlib.pyplot as plt
import odc.stac
from shapely.geometry import Point,Polygon
import xrspatial

catalog = pystac_client.Client.open(
    "https://planetarycomputer.microsoft.com/api/stac/v1",
    modifier=planetary_computer.sign_inplace,
)


def add_data(df, output_cols, collection, bands, time_range,crs, derived_bands = None,odc_args = dict(), 
             verbose = 2, shapely_in_bds = False):
    start = time.time()
    for output_col in output_cols:
        df[output_col] = np.nan
    visited = pd.Series(False, index = range(len(df)))
    
    if derived_bands is not None:
        all_bands = bands + derived_bands
    else:
        all_bands = bands
        
    while(sum(~visited) > 0):
        asset_start = time.time()
        first_unvisited_idx = np.min(np.where(~visited))
        if verbose >= 2:
            print("Pulling asset for sample point at index ", first_unvisited_idx, "...")
        target_pt = df.iloc[first_unvisited_idx]
        point = {"type": "Point", "coordinates": [target_pt['lon'], target_pt['lat']]}
        item, asset = fetch_asset(collection, bands, point, time_range, crs, odc_args)
        
        if verbose >=2:
            print(f"Asset pulled in {time.time() - asset_start} seconds")
        
        if derived_bands is not None:
            for band in derived_bands:
                asset = derived_bands_fns[band](asset)
        
        if shapely_in_bds:
            poly = Polygon(item.geometry["coordinates"][0])
            is_in_bds = lambda x,y: poly.contains(Point(x,y))
            in_bds_idx = [is_in_bds(df['lon'][i],df['lat'][i]) for i in range(len(df))]
        else:
            is_in_bds = get_is_in_bds(asset)
            in_bds_idx = is_in_bds(df['lon'],df['lat'])
        in_bds = df[in_bds_idx]
        visited[in_bds_idx] = True
        
        if verbose >=2:
            print(str(len(in_bds)), "points within new asset")
    
        for i,band in enumerate(all_bands):
            pts = zip(in_bds["lon"], in_bds["lat"])
            data = []
            for p in pts:
                data.append(asset.sel({"lon":p[0],"lat":p[1]},method = "nearest")[band].item())
            df.loc[np.where(in_bds_idx)[0],output_cols[i]] = data
        
        if verbose >= 2:
            print(f'{sum(visited)}/{len(visited)}')
            print(f"Asset processed in {time.time() - asset_start} seconds")

    if verbose >=1:
        print(f"Sampling completed in {time.time() - start} seconds")
    
    return(df)

def add_aspect(asset, dem_band='data'):
    aspects = xrspatial.aspect(asset[dem_band].squeeze())[1:-1,1:-1]
    aspects.pad(lon = (1,1),lat = (1,1))
    asset['aspect'] = aspects
    return(asset)

def add_grade(asset, dem_band='data'):
    grad_x,grad_y = np.gradient(asset[dem_band].squeeze())
    grad = np.sqrt(grad_x**2 + grad_y**2) / 30
    asset['grade_lon'] = (['lat','lon'],grad_x)
    asset['grade_lat'] = (['lat','lon'],grad_y)
    asset['grade'] = (['lat','lon'],grad)
    asset['grade'] = asset['grade'].rolling(lat = 3,lon = 3, 
                                        center = True, min_periods = 1).mean()
    return(asset)

derived_bands_fns = {"aspect":add_aspect,
                     "grade":add_grade}

def add_ndvis(sample_pts, min_idx, roll_window, crs):

    null_idxs = sample_pts.index[sample_pts['ndvi'].isnull()]
    null_idxs = null_idxs[null_idxs > min_idx]
    target_pt_idx = np.min(null_idxs)
    print("Pulling asset for sample point at index ", target_pt_idx, "...")
    
    target_pt = sample_pts.iloc[target_pt_idx]
    point = {"type": "Point", "coordinates": [target_pt['lon'], target_pt['lat']]}
    data = fetch_asset("landsat-c2-l2", point, crs)
    data = calc_ndvi(data, roll_window)

    is_in_bds = get_is_in_bds(data)
    in_bds_idx = is_in_bds(sample_pts['lon'],sample_pts['lat'])
    in_bds = sample_pts[in_bds_idx]
    print(str(len(in_bds)), "points within new asset")
    
    pts = zip(in_bds["lon"], in_bds["lat"])
    ndvis = []
    for p in pts:
        ndvis.append(data.sel({"lon":p[0],"lat":p[1]},method = "nearest")["log_ndvi"].item())
    sample_pts.loc[np.where(in_bds_idx)[0],'ndvi'] = ndvis
    
    non_nan_num = (~sample_pts['ndvi'].isnull()).sum()
    total_len = len(sample_pts['ndvi'])
    n_left = total_len - non_nan_num
    print(f'{non_nan_num}/{total_len}')

    return(target_pt_idx, sample_pts,n_left)

def calc_ndvi(data, roll_window):
    nir_arr = data["nir08"]
    red_arr = data["red"]
    denom_arr = nir_arr + red_arr
    data["ndvi"] = (nir_arr - red_arr) / (nir_arr + red_arr)
    data["ndvi"] = data["ndvi"].where(denom_arr != 0, 0)
    data["ndvi"] = data["ndvi"].where(data["ndvi"]<10,10)
    data["log_ndvi"] = np.log(data["ndvi"]+1)
    data["ndvi"] = data["ndvi"].rolling(lat = roll_window,lon = roll_window, 
                                        center = True, min_periods = 1).mean()
    data["log_ndvi"] = data["log_ndvi"].rolling(lat = roll_window,lon = roll_window, 
                                            center = True, min_periods = 1).mean()
    return(data)

def fetch_asset(collection, bands, target_pt, time_range, crs, odc_args = dict()):
    search = catalog.search(collections=[collection], intersects=target_pt, datetime=time_range)
    items = search.item_collection()
    data = odc.stac.load(
        items[:1], bands=bands,intersects=target_pt, **odc_args
    )#.squeeze()
    data = data.rio.reproject(crs)
    data = data.rename({"x":"lon","y":"lat"})
    return(items[0], data)

def get_bounds(sum_image):
    nz_lon_idx = sum_image.sum(dim = 'lat').to_numpy().nonzero()
    nz_lat_idx = sum_image.sum(dim = 'lon').to_numpy().nonzero()

    mnlnidx = np.min(nz_lon_idx); mxlnidx = np.max(nz_lon_idx);
    mnltidx = np.min(nz_lat_idx); mxltidx = np.max(nz_lat_idx);

    lb_rb = sum_image['lon'][[mnlnidx,mxlnidx]]
    bb_tb = sum_image['lat'][[mnltidx, mxltidx]]
    left_bound = min(lb_rb).item(); right_bound = max(lb_rb).item()
    bottom_bound = min(bb_tb).item(); top_bound = max(bb_tb).item()
    
    return((left_bound, right_bound, top_bound, bottom_bound))

def get_is_in_bds(data_xr, return_lines = False):
    
    sum_image = get_sum_image(data_xr)
    left_vertex, right_vertex, top_vertex, bottom_vertex = get_vertices(sum_image)
    
    #m = (y2-y1)/(x2-x1)
    lv_tv_slope = (top_vertex[1] - left_vertex[1])/(top_vertex[0] - left_vertex[0])
    lv_bv_slope = (bottom_vertex[1] - left_vertex[1])/(bottom_vertex[0] - left_vertex[0])
    
    #y=mx+b -> b = y - mx
    lv_tv_yintercept = left_vertex[1] - (lv_tv_slope * left_vertex[0])
    bv_rv_yintercept = bottom_vertex[1] - (lv_tv_slope * bottom_vertex[0])
    lv_bv_yintercept = left_vertex[1] - (lv_bv_slope * left_vertex[0])
    tv_rv_yintercept = top_vertex[1] - (lv_bv_slope * top_vertex[0])

    def is_in_bds(lon, lat):
        lv_tv_transform = lat - (lv_tv_slope * lon)
        lv_bv_transform = lat - (lv_bv_slope * lon)
        in_ud_dir = (lv_tv_transform <= lv_tv_yintercept) & (lv_tv_transform >= bv_rv_yintercept)
        in_lr_dir = (lv_bv_transform <= tv_rv_yintercept) & (lv_bv_transform >= lv_bv_yintercept)
        return(in_ud_dir & in_lr_dir)
    
    if not return_lines:
        return(is_in_bds)
    else:
        lines = {"m1":lv_tv_slope,"b1":lv_tv_yintercept,
                 "m2":lv_tv_slope,"b2":bv_rv_yintercept,
                 "m3":lv_bv_slope,"b3":lv_bv_yintercept,
                 "m4":lv_bv_slope,"b4":tv_rv_yintercept}
        return(is_in_bds, lines)
    
def get_sum_image(data_xr):
    images = [data_xr[band] for band in data_xr]
    sum_image = images[0]
    for img in images[1:]:
        sum_image = sum_image + img

    return(sum_image)

def get_vertices(sum_image):
    
    lb, rb, ub, bb = get_bounds(sum_image)
    
    left_col_nz = sum_image.loc[:,lb].to_numpy().nonzero()[0][0]
    right_col_nz = sum_image.loc[:,rb].to_numpy().nonzero()[0][0]
    top_row_nz = sum_image.loc[ub,:].to_numpy().nonzero()[0][0]
    bottom_row_nz = sum_image.loc[bb,:].to_numpy().nonzero()[0][0]

    left_vertex_lat = sum_image['lat'][left_col_nz].item()
    right_vertex_lat = sum_image['lat'][right_col_nz].item()
    top_vertex_lon = sum_image['lon'][top_row_nz].item()
    bottom_vertex_lon = sum_image['lon'][bottom_row_nz].item()

    left_vertex = (lb, left_vertex_lat)
    right_vertex = (rb, right_vertex_lat)
    top_vertex = (top_vertex_lon, ub)
    bottom_vertex = (bottom_vertex_lon, bb)
    
    return((left_vertex, right_vertex, top_vertex, bottom_vertex))
    