# -*- coding: utf-8 -*-

import overpy

api = overpy.Overpass()

# fetch all ways and nodes
result = api.query("""
    way(50.746,7.154,50.748,7.157) ["highway"];
    (._;>;);
    out body;
    """)

for way in result.ways:
    print("Name: %s" % way.tags.get("name", "n/a"))
    print("  Highway: %s" % way.tags.get("highway", "n/a"))
    print("  Nodes:")
    for node in way.nodes:
        print("    Lat: %f, Lon: %f" % (node.lat, node.lon))
        
way(poly:"32.1085856 34.8051541 32.1087969 34.8052841 32.1086572 34.8056320 32.1090022 34.8052390 32.1086452 34.8050192")["building" = "yes"];
/*added by auto repair*/
(._;>;);
/*end of auto repair*/
out;