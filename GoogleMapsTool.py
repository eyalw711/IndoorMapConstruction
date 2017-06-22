import gmplot


class GoogleMapsConfig:
    def __init__(self, latPoints, longPoints, fileName):
        self.Lat = latPoints
        self.Long = longPoints
        self.FileName = fileName


class GoogleMapFactory:
    @staticmethod
    def GenerateMap(configuration, marker=False, size=0.5):
        assert isinstance(configuration, GoogleMapsConfig)
        print("Starting generate Google maps plot to file {}".format(configuration.FileName))
        assert len(configuration.Lat) == len(configuration.Long)
        if len(configuration.Lat) <= 0:
            print("Cannot create google map with empty lat/long list")
            return
        gomap = gmplot.GoogleMapPlotter(configuration.Lat[0], configuration.Long[0], 20)
        gomap.scatter(configuration.Lat, configuration.Long, marker=marker, size=size)
        gomap.draw("GoogleMaps//{}".format(configuration.FileName))