'''
Modified from linqu's answer on stackoverflow question:
how-to-convert-from-longitude-and-latitude-to-country-or-city
'''
import requests
from shapely.geometry import mapping, shape
from shapely.prepared import prep
from shapely.geometry import Point

data = requests.get("https://raw.githubusercontent.com/datasets/geo-countries/master/data/countries.geojson").json()

countries = {}
for feature in data["features"]:
    geom = feature["geometry"]
    country = feature["properties"]["ADMIN"]
    countries[country] = prep(shape(geom))

def get_country(lon, lat):
    ''' get country name from long lat
    '''
    point = Point(lon, lat)
    for country, geom in countries.items():
        if geom.contains(point):
            return country

    return "NA"
