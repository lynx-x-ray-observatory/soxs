import numpy as np

class SpatialModel(object):
    def __init__(self, ra, dec):
        self.ra = ra
        self.dec = dec
        self.num_events = self.ra.size

    def __add__(self, other):
        ra = np.concatenate([self.ra, other.ra])
        dec = np.concatenate([self.dec, other.dec])
        return SpatialModel(ra, dec)

class PointSourceModel(SpatialModel):
    def __init__(self, pt_ra, pt_dec, num_events):
        ra = pt_ra*np.ones(num_events)
        dec = pt_dec*np.ones(num_events)
        super(PointSourceModel, self).__init__(ra, dec)
