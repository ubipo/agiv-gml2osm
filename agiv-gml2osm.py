#!/usr/bin/env python3

from math import *
import xml.etree.ElementTree as et
import argparse
import os
import progressbar

parser = argparse.ArgumentParser(description='Convert AGIV GRB gml XML files with lambert72 coordinates to the osm XML format.')
parser.add_argument('src', help='the directory to scan or .gml file to convert')
parser.add_argument('des', nargs='?', default=os.getcwd(), help='the .osm file to write to')
parser.add_argument('-v', '--verbose', help='run in verbose (debug) mode')
args = parser.parse_args()

def lambert72toWGS84(x, y):
    n = 0.77164219
    F = 1.81329763
    thetaFudge = 0.00014204
    e = 0.08199189
    a = 6378388
    xDiff = 149910
    yDiff = 5400150
    theta0 = 0.07604294

    xReal = xDiff - x
    yReal = yDiff - y

    rho = sqrt(xReal * xReal + yReal * yReal)
    theta = atan(xReal / -yReal)

    longitude = (theta0 + (theta + thetaFudge) / n) * 180 / pi;
    latitude = 0;

    for _ in range(5):
        latitude = (2 * atan(pow(F * a / rho, 1 / n) * pow((1 + e * sin(latitude)) / (1 - e * sin(latitude)), e / 2))) - pi / 2;

    latitude *= 180 / pi;

    return (latitude, longitude);


def haversine(lat1, lon1, lat2, lon2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    r = 6371008 # Radius of earth in metres.
    return c * r


def findCloseNode(coords, nodes):
    """ 
    Check if a node is already in the vicinity, if yes reuse.
    """
    for node in nodes:
        # Use haversine to calculate distance
        if haversine(*coords, float(node['attrib']['lat']), float(node['attrib']['lon'])) < 0.2:
            return node
    return None


# if args.src[-4:] == '.gml':

# Parse
gml = et.parse(args.src)
root = gml.getroot()

nodes = []
ways = []
relations = []
osmId = -70000

features = root.findall('.//{http://www.opengis.net/gml}featureMember')

with progressbar.ProgressBar(max_value=len(features)) as bar:
    for i, feature in enumerate(features):
        featureType = feature.find('.//{http://www.agiv.be/agiv}TYPE').text
        featureOsmType = ''
        if featureType == '1':
            featureOsmType = 'house'
        elif featureType == '2':
            featureOsmType = 'detached'
        polygon = feature.find('.//{http://www.opengis.net/gml}Polygon')
        outer = polygon.findall('{http://www.opengis.net/gml}outerBoundaryIs')
        inner = polygon.findall('{http://www.opengis.net/gml}innerBoundaryIs')
        for boundary in outer+inner:
            # Scan for and split node coordinates
            coordTexts = boundary.find('.//{http://www.opengis.net/gml}coordinates').text.split(' ')

            # Split in lat and lon floats
            coordsLambert72 = [(float(comp) for comp in coordText.split(',')) for coordText in coordTexts]

            # Convert whole list to WGS84
            coords = [lambert72toWGS84(*coordLambert72) for coordLambert72 in coordsLambert72]

            # Simple checks
            if coords[0] != coords[-1]:
                print('Unclosed way! Skipping...')
                break

            if len(coords) <= 2:
                print('One long way! Skipping...')
                break

            # Give each osm elem a unique negative id
            osmId -= 1

            # Create a new way
            way = {
                'attrib': {
                    'id': str(osmId),
                    'action': 'modify'
                },
                'nd': [],
                'tag': [{
                    'attrib': {
                        'k': 'building',
                        'v': featureOsmType
                    }
                }]
            }

            # Create all the nodes (one for each coordinate pair), and add them to the way
            for j, coord in enumerate(coords):
                closeNode = findCloseNode(coord, nodes)
                if closeNode:
                    nodeId = closeNode['attrib']['id']
                    # A node can't follow itself up
                    if nodeId in lastNodeId:
                        continue
                else:
                    osmId -= 1
                    nodeId = str(osmId)
                    lastNodeId = nodeId
                    nodes.append({
                        'attrib': {
                            'id': str(osmId),
                            'action': 'modify',
                            'lat': str(coord[0]),
                            'lon': str(coord[1])
                        }
                    })

                way['nd'].append({
                    'attrib': {
                        'ref': nodeId,
                    }
                })
                

            osmId -= 1
            ways.append(way)
        # / inner + outer
        bar.update(i)
    # / features


osm = et.ElementTree()
osm._setroot(et.Element('osm', attrib = {
        'version':'0.6',
        'generator':'agiv-gml2osm'
    }))
osm_root = osm.getroot()

for node in nodes:
    osm_root.append(et.Element('node', attrib = node['attrib']))

for way in ways:
    osm_root.append(et.Element('way', attrib = way['attrib']))
    way_xml = osm_root.find(".//*[@id='"+way['attrib']['id']+"']")
    for nd in way['nd']:
        way_xml.append(et.Element('nd', attrib = nd['attrib']))
    for tag in way['tag']:
        way_xml.append(et.Element('tag', attrib = tag['attrib']))

osm.write('output.osm')
