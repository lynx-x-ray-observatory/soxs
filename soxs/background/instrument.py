from soxs.utils import parse_prng, \
    parse_value, mylog, create_region, get_data_file
from soxs.background.events import make_diffuse_background
import numpy as np
import astropy.io.fits as pyfits
from regions import PixCoord


class InstrumentalBackground:

    def __init__(self, channel, count_rate, default_focal_length,
                 exp_time):
        self.channel = channel
        self.count_rate = count_rate
        self.default_focal_length = default_focal_length
        self.exp_time = exp_time

    @classmethod
    def from_filename(cls, filename, ext_area, focal_length):
        """
        Read an instrumental background spectrum from 
        a FITS PHA file. 

        Parameters
        ----------
        filename : string
            The path to the file containing the spectrum.
        focal_length : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The default focal length of the instrument
            in meters. 
        """
        fn = get_data_file(filename)
        with pyfits.open(fn) as f:
            hdu = f["SPECTRUM"]
            exp_time = hdu.header["EXPOSURE"]
            if "COUNTS" in hdu.data.names:
                count_rate = hdu.data["COUNTS"]/exp_time
            else:
                count_rate = hdu.data["COUNT_RATE"]
            count_rate /= ext_area
            channel = hdu.data["CHANNEL"]
        return cls(channel, count_rate, focal_length,
                   exp_time)

    def generate_channel_spectrum(self, t_exp, solid_angle, 
                                  focal_length=None, prng=None):
        """
        Generate photon energy channels from this instrumental
        background spectrum given an exposure time,
        effective area, and solid angle.

        Parameters
        ----------
        t_exp : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The exposure time in seconds.
        solid_angle : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The solid angle in arcmin**2.
        focal_length : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
            The focal length in meters. Default is to use
            the default focal length of the instrument
            configuration.
        prng : :class:`~numpy.random.RandomState` object, integer, or None
            A pseudo-random number generator. Typically will only 
            be specified if you have a reason to generate the same 
            set of random numbers, such as for a test. Default is None, 
            which sets the seed based on the system time. 
        """
        t_exp = parse_value(t_exp, "s")
        solid_angle = parse_value(solid_angle, "arcmin**2")
        prng = parse_prng(prng)
        if focal_length is None:
            focal_length = self.default_focal_length
        else:
            focal_length = parse_value(focal_length, "m")
        fac = t_exp * solid_angle # Backgrounds are normalized to 1 arcmin**2
        fac *= (focal_length / self.default_focal_length) ** 2
        return prng.poisson(lam=self.count_rate*fac)

    def generate_channels(self, t_exp, solid_angle, focal_length=None,
                          prng=None):
        """
        Generate photon energy channels from this instrumental
        background spectrum given an exposure time,
        effective area, and solid angle.

        Parameters
        ----------
        t_exp : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The exposure time in seconds.
        solid_angle : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The solid angle in arcmin**2.
        focal_length : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
            The focal length in meters. Default is to use
            the default focal length of the instrument
            configuration.
        prng : :class:`~numpy.random.RandomState` object, integer, or None
            A pseudo-random number generator. Typically will only 
            be specified if you have a reason to generate the same 
            set of random numbers, such as for a test. Default is None, 
            which sets the seed based on the system time. 
        """
        ncts = self.generate_channel_spectrum(t_exp, solid_angle,
                                              focal_length=focal_length, 
                                              prng=prng)
        return np.concatenate([self.channel[i]*np.ones(n)
                               for i, n in enumerate(ncts)])


def make_instrument_background(inst_spec, event_params, rmf, prng=None):
    from collections import defaultdict

    prng = parse_prng(prng)

    bkgnd_spec = inst_spec["bkgnd"]

    if isinstance(bkgnd_spec[0], str):
        nchips = len(event_params["chips"])
        bkgnd_spec = [bkgnd_spec]*nchips

    bkg_events = defaultdict(list)
    pixel_area = (event_params["plate_scale"]*60.0)**2
    for i, chip in enumerate(event_params["chips"]):
        rtype = chip[0]
        args = chip[1:]
        r, bounds = create_region(rtype, args, 0.0, 0.0)
        sa = (bounds[1]-bounds[0])*(bounds[3]-bounds[2])*pixel_area
        bspec = InstrumentalBackground.from_filename(
            bkgnd_spec[i][0], bkgnd_spec[i][1], inst_spec['focal_length'])
        chan = bspec.generate_channels(
            event_params["exposure_time"], sa, prng=prng)
        n_events = chan.size
        detx = prng.uniform(low=bounds[0], high=bounds[1], size=n_events)
        dety = prng.uniform(low=bounds[2], high=bounds[3], size=n_events)
        if rtype in ["Box", "Rectangle"]:
            thisc = slice(None, None, None)
            n_det = n_events
        else:
            thisc = r.contains(PixCoord(detx, dety))
            n_det = thisc.sum()
        ch = chan[thisc].astype('int')
        e = rmf.ch_to_eb(ch, prng=prng)
        bkg_events["energy"].append(e)
        bkg_events[rmf.chan_type].append(ch)
        bkg_events["detx"].append(detx[thisc])
        bkg_events["dety"].append(dety[thisc])
        bkg_events["chip_id"].append(i*np.ones(n_det))
    for key in bkg_events:
        bkg_events[key] = np.concatenate(bkg_events[key])

    if bkg_events["energy"].size == 0:
        raise RuntimeError("No instrumental background events were detected!!!")
    else:
        mylog.info(f"Making {bkg_events['energy'].size} events "
                   f"from the instrumental background.")

    return make_diffuse_background(bkg_events, event_params, rmf, prng=prng)
