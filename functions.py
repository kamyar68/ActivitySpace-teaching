import geopandas as gpd
from shapely.geometry import LineString, Point

def generate_points_along_routes(input_gdf, output_path, distance_interval=30):
    """
    Generate points along routes in a GeoDataFrame at specified distance intervals.

    Parameters:
    - input_gdf (geopandas.GeoDataFrame): Input GeoDataFrame containing route geometries.
    - output_path (str): Path to save the resulting GeoDataFrame with points.
    - distance_interval (float): Distance interval between generated points along each route.

    Returns:
    - geopandas.GeoDataFrame: GeoDataFrame containing points along routes.
    """
    # Create an empty list to store rows
    points_rows = []

    for index, row in input_gdf.iterrows():
        # Get the route geometry
        route_geometry = row['geometry']

        # Generate points along the route at the specified interval
        points_along_route = []
        length = route_geometry.length

        # Ensure the route has sufficient length for generating points
        if length >= distance_interval:
            current_distance = 0

            while current_distance < length:
                point_on_route = route_geometry.interpolate(current_distance)

                # Create a new row for each point
                new_row = row.copy()
                new_row['geometry'] = point_on_route
                points_rows.append(new_row)

                current_distance += distance_interval

    # Create the GeoDataFrame from the list of rows
    points_gdf = gpd.GeoDataFrame(points_rows, geometry='geometry', crs=input_gdf.crs)

    # Save the resulting GeoDataFrame to a file
    #points_gdf.to_file(output_path)

    return points_gdf


def generate_points_along_polygon(input_gdf, output_path, distance_interval=30):
    """
    Generate points along the exterior outline of Polygon geometries in a GeoDataFrame at specified distance intervals.

    Parameters:
    - input_gdf (geopandas.GeoDataFrame): Input GeoDataFrame containing Polygon geometries.
    - output_path (str): Path to save the resulting GeoDataFrame with points.
    - distance_interval (float): Distance interval between generated points along each polygon outline.

    Returns:
    - geopandas.GeoDataFrame: GeoDataFrame containing points along polygon outlines.
    """
    # Create an empty list to store rows
    points_rows = []

    for index, row in input_gdf.iterrows():
        # Get the polygon geometry
        polygon = row['geometry']

        # Generate points along the exterior outline of the polygon at the specified interval
        points_along_outline = []
        length = polygon.exterior.length

        # Ensure the polygon has sufficient length for generating points
        if length >= distance_interval:
            current_distance = 0

            while current_distance < length:
                point_on_outline = polygon.exterior.interpolate(current_distance)

                # Create a new row for each point
                new_row = row.copy()
                new_row['geometry'] = Point(point_on_outline.x, point_on_outline.y)
                points_rows.append(new_row)

                current_distance += distance_interval

    # Create the GeoDataFrame from the list of rows
    points_gdf = gpd.GeoDataFrame(points_rows, geometry='geometry', crs=input_gdf.crs)

    # Save the resulting GeoDataFrame to a file
    #points_gdf.to_file(output_path)

    return points_gdf

