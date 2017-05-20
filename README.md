# Indoor Map Construction From Crowd-Source  
This is an final project for a B.Sc degree in Electrical Engineering in Tel-Aviv University.  
In our project the Map-Constructor will build navigable indoor map objects of the buildings floors  
given a data-base of pedestrians trajectories.

Check the [Wiki](https://github.com/eyalw711/IndoorMapConstruction/wiki) for more details.

### Last Results:  
Clustering simulated trajectories in a "W" shaped building:  
![alt text](https://github.com/eyalw711/IndoorMapConstruction/blob/master/results/Mapping%20Simulated%20DB%20with%20eps%201-2%20minlns%205.png "Clustering Result")
![alt text](https://github.com/eyalw711/IndoorMapConstruction/blob/master/results/basicMap.png "Basic Map")


### Installation Instructions:  
Work with anaconda environment with python 3.6

Important packages:

- geopy  
A Geo-Navigation Package  
[Documentation](https://pypi.python.org/pypi/geopy)  
To install: `pip install geopy`

- nvector  
A Geo-Navigation Package  
[Documentation](https://pypi.python.org/pypi/nvector)  
To install: `pip install nvector`  

- shapely  
A Geometry Package  
[Documentation](http://toblerity.org/shapely/manual.html)  
To install: `conda install -c conda-forge shapely=1.5.17`  

- overpy  
An Open Street Map Interface Package  
[Documentation](https://pypi.python.org/pypi/overpy/)  
To install: `pip install overpy`

- networkx  
Graphs Algorithms Package  
[Documentation](http://networkx.readthedocs.io/en/networkx-1.11/index.html)  
To install: `pip install networkx`

- pickle  
Python Objects Serialization  
[Documentation](https://docs.python.org/3/library/pickle.html)  
To install: `pip install pickle`
