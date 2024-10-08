import numpy as np
from astropy.io import fits
from pathlib import PurePath

from soxs.events.utils import _region_filter
from soxs.utils import parse_value


def _write_spectrum(
    bins,
    spec,
    parameters,
    specfile,
    overwrite=False,
    noisy=True,
):
    exp_time = parameters.get("EXPOSURE", None)
    spectype = parameters["CHANTYPE"]

    if noisy:
        cnt_fmt = "1J"
        cnt_type = "int32"
    else:
        cnt_fmt = "1D"
        cnt_type = "float64"

    col_ch = fits.Column(name="CHANNEL", format="1J", array=bins)
    col_chf = fits.Column(
        name=spectype.upper(), format="1D", array=bins.astype("float64")
    )

    cols = [col_ch, col_chf]

    if exp_time is None:
        rate = spec
    else:
        rate = spec / exp_time
        col_cnt = fits.Column(
            name="COUNTS", format=cnt_fmt, array=spec.astype(cnt_type)
        )
        cols.append(col_cnt)

    col_rate = fits.Column(name="COUNT_RATE", format="1D", array=rate)
    cols.append(col_rate)

    coldefs = fits.ColDefs(cols)

    tbhdu = fits.BinTableHDU.from_columns(coldefs)
    tbhdu.name = "SPECTRUM"

    tbhdu.header["DETCHANS"] = spec.size
    tbhdu.header["TOTCTS"] = spec.sum()
    tbhdu.header["EXPOSURE"] = exp_time
    tbhdu.header["LIVETIME"] = exp_time
    tbhdu.header["CONTENT"] = spectype
    tbhdu.header["HDUCLASS"] = "OGIP"
    tbhdu.header["HDUCLAS1"] = "SPECTRUM"
    tbhdu.header["HDUCLAS2"] = "TOTAL"
    tbhdu.header["HDUCLAS3"] = "TYPE:I"
    tbhdu.header["HDUCLAS4"] = "COUNT"
    tbhdu.header["HDUVERS"] = "1.1.0"
    tbhdu.header["HDUVERS1"] = "1.1.0"
    tbhdu.header["CHANTYPE"] = spectype
    tbhdu.header["BACKFILE"] = "none"
    tbhdu.header["CORRFILE"] = "none"
    tbhdu.header["POISSERR"] = noisy
    for key in ["RESPFILE", "ANCRFILE", "MISSION", "TELESCOP", "INSTRUME"]:
        tbhdu.header[key] = parameters[key]
    tbhdu.header["AREASCAL"] = 1.0
    tbhdu.header["CORRSCAL"] = 0.0
    tbhdu.header["BACKSCAL"] = 1.0

    hdulist = fits.HDUList([fits.PrimaryHDU(), tbhdu])

    hdulist.writeto(specfile, overwrite=overwrite)


def _make_spectrum(
    evtfile,
    region=None,
    format="ds9",
    exclude=False,
    emin=None,
    emax=None,
    tmin=None,
    tmax=None,
):
    from soxs.response import RedistributionMatrixFile

    parameters = {}
    if isinstance(evtfile, str):
        with fits.open(evtfile) as f:
            hdu = f["EVENTS"]
            evt_mask = np.ones(hdu.data["ENERGY"].size, dtype="bool")
            if region is not None:
                evt_mask &= _region_filter(hdu, region, format=format, exclude=exclude)
            if tmin is not None:
                tmin = parse_value(tmin, "s")
                evt_mask &= hdu.data["TIME"] > tmin
            else:
                tmin = 0.0
            if tmax is not None:
                tmax = parse_value(tmax, "s")
                evt_mask &= hdu.data["TIME"] < tmax
            else:
                tmax = hdu.header["EXPOSURE"]
            if emin is not None:
                emin = parse_value(emin, "keV") * 1000.0
                evt_mask &= hdu.data["ENERGY"] > emin
            if emax is not None:
                emax = parse_value(emax, "keV") * 1000.0
                evt_mask &= hdu.data["ENERGY"] < emax
            spectype = hdu.header["CHANTYPE"]
            rmf = hdu.header["RESPFILE"]
            p = hdu.data[spectype][evt_mask]
            parameters["EXPOSURE"] = tmax - tmin
            for key in ["RESPFILE", "ANCRFILE", "MISSION", "TELESCOP", "INSTRUME"]:
                parameters[key] = hdu.header[key]
    else:
        rmf = evtfile["rmf"]
        spectype = evtfile["channel_type"]
        p = evtfile[spectype]
        parameters["RESPFILE"] = PurePath(rmf).parts[-1]
        parameters["ANCRFILE"] = PurePath(evtfile["arf"]).parts[-1]
        parameters["TELESCOP"] = evtfile["telescope"]
        parameters["INSTRUME"] = evtfile["instrument"]
        parameters["MISSION"] = evtfile["mission"]
        parameters["EXPOSURE"] = evtfile["exposure_time"]
    parameters["CHANTYPE"] = spectype

    rmf = RedistributionMatrixFile(rmf)
    spec = np.bincount(p, minlength=rmf.n_ch + rmf.cmin)[rmf.cmin :]
    bins = (np.arange(rmf.n_ch) + rmf.cmin).astype("int32")

    return bins, spec, parameters


def write_spectrum(
    evtfile,
    specfile,
    region=None,
    format="ds9",
    exclude=False,
    emin=None,
    emax=None,
    tmin=None,
    tmax=None,
    overwrite=False,
):
    r"""
    Bin event energies into a spectrum and write it to
    a FITS binary table. Does not do any grouping of
    channels, and will automatically determine PI or PHA.
    Optionally allows filtering based on time, energy, or
    spatial region.

    Parameters
    ----------
    evtfile : string
        The name of the event file to read the events from.
    specfile : string
        The name of the spectrum file to be written.
    region : string, Region, or Regions, optional
        The region(s) to be used for the filtering. Default: None
    format : string, optional
        The file format specifier for the region. "ds9",
        "crtf", "fits", etc. Default: "ds9"
    exclude : boolean, optional
        If True, the events in a specified *region* will be excluded
        instead of included in the spectrum. Default: False
    emin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The minimum energy of the events to be included, in keV.
        Default is the lowest energy available.
    emax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The maximum energy of the events to be included, in keV.
        Default is the highest energy available.
    tmin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The minimum energy of the events to be included, in seconds.
        Default is the earliest time available.
    tmax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The maximum energy of the events to be included, in seconds.
        Default is the latest time available.
    overwrite : boolean, optional
        Whether to overwrite an existing file with
        the same name. Default: False
    """

    bins, spec, parameters = _make_spectrum(
        evtfile,
        region=region,
        format=format,
        exclude=exclude,
        emin=emin,
        emax=emax,
        tmin=tmin,
        tmax=tmax,
    )
    _write_spectrum(bins, spec, parameters, specfile, overwrite=overwrite)


def plot_spectrum(
    specfile,
    plot_energy=True,
    ebins=None,
    lw=2,
    xmin=None,
    xmax=None,
    ymin=None,
    ymax=None,
    xscale=None,
    yscale=None,
    label=None,
    fontsize=18,
    fig=None,
    ax=None,
    plot_counts=False,
    per_keV=True,
    noerr=False,
    plot_used=False,
    **kwargs,
):
    """
    Make a quick Matplotlib plot of a convolved spectrum
    from a file. A Matplotlib figure and axis is returned.

    Parameters
    ----------
    specfile : string
        The file to be opened for plotting.
    plot_energy : boolean, optional
        Whether to plot in energy or channel space. Default is
        to plot in energy, unless the RMF for the spectrum
        cannot be found.
    ebins : integer, tuple, or NumPy array, optional
        If set, these are the energy bin edges in which the spectrum
        will be binned. If an integer, the counts spectrum will be reblocked
        by this number. If a 2-tuple, the first element is the minimum
        significance (assuming Poisson statistics) of each bin and the
        second element is the maximum number of channels to be combined
        in the bin. If a NumPy array, these are the bins that will be
        used. If not set, the counts will be binned according
        to channel. Default: None
    lw : float, optional
        The width of the lines in the plots. Default: 2.0 px.
    xmin : float, optional
        The left-most energy (in keV) or channel to plot. Default is the
        minimum value in the spectrum.
    xmax : float, optional
        The right-most energy (in keV) or channel to plot. Default is the
        maximum value in the spectrum.
    ymin : float, optional
        The lower extent of the y-axis. By default it is set automatically.
    ymax : float, optional
        The upper extent of the y-axis. By default it is set automatically.
    xscale : string, optional
        The scaling of the x-axis of the plot. Default: "linear"
    yscale : string, optional
        The scaling of the y-axis of the plot. Default: "linear"
    label : string, optional
        The label of the spectrum. Default: None
    fontsize : int
        Font size for labels and axes. Default: 18
    fig : :class:`~matplotlib.figure.Figure`, optional
        A Figure instance to plot in. Default: None, one will be
        created if not provided.
    ax : :class:`~matplotlib.axes.Axes`, optional
        An Axes instance to plot in. Default: None, one will be
        created if not provided.
    plot_counts : boolean, optional
        If set to True, the counts instead of the count rate will
        be plotted. Default: False
    per_keV : boolean, optional
        If set to True, the spectrum will be plotted in counts per keV.
        Otherwise, it will be plotted in counts per bin. Default: True
    noerr : boolean, optional
        If True, the spectrum will be plotted without errorbars. This
        will always be the case if the spectrum is not noisy.
        Default: False
    plot_used : boolean, optional
        If set to True, only the bins which contain more than 0
        counts will be plotted. Default: False

    Returns
    -------
    A tuple of the :class:`~matplotlib.figure.Figure` object, the
    :class:`~matplotlib.axes.Axes` object, and the NumPy array of
    energy bins that are used.
    """
    import matplotlib.pyplot as plt
    from astropy.nddata import block_reduce

    from soxs.instrument import RedistributionMatrixFile

    f = fits.open(specfile)
    hdu = f["SPECTRUM"]
    chantype = hdu.header["CHANTYPE"]
    y = hdu.data["COUNTS"].astype("float64")
    if not hdu.header["POISSERR"]:
        noerr = True
    if plot_energy:
        rmf = hdu.header.get("RESPFILE", None)
        if rmf is not None:
            rmf = RedistributionMatrixFile(rmf)
            emin = rmf.ebounds_data["E_MIN"]
            emax = rmf.ebounds_data["E_MAX"]
            e = 0.5 * (emin + emax)
            if ebins is None:
                ebins = np.append(emin, emax[-1])
            elif isinstance(ebins, int):
                y = block_reduce(y, ebins)
                ebins = np.append(emin[::ebins], emax[-1])
            else:
                if isinstance(ebins, tuple):
                    if len(ebins) != 2:
                        raise RuntimeError(
                            "If ebins is a tuple, it must have 2 elements!"
                        )
                    sigma, max_size = ebins
                    sigma2 = sigma * sigma
                    ebins = [emin[0]]
                    sum = 0.0
                    bin_size = 0
                    for i in range(y.size):
                        sum += y[i]
                        bin_size += 1
                        if sum >= sigma2 or bin_size == max_size or i == y.size - 1:
                            ebins.append(emax[i])
                            sum = 0.0
                            max_size = 0
                    ebins = np.array(ebins)
                y = np.histogram(e, ebins, weights=y)[0].astype("float64")
            xlabel = "Energy (keV)"
            xmid = 0.5 * (ebins[1:] + ebins[:-1])
            xerr = 0.5 * np.diff(ebins)
        else:
            raise RuntimeError(
                "Cannot find the RMF associated with this "
                "spectrum, so I cannot plot in energy!"
            )
    else:
        xmid = hdu.data[chantype]
        xerr = 0.5
        xlabel = f"Channel ({chantype})"
    dx = 2.0 * xerr
    yerr = np.sqrt(y)
    if not plot_counts:
        y /= hdu.header["EXPOSURE"]
        yerr /= hdu.header["EXPOSURE"]
    if per_keV and plot_energy:
        yunit = "keV"
        y /= dx
        yerr /= dx
    else:
        yunit = "bin"
    f.close()
    if fig is None:
        fig = plt.figure(figsize=(10, 10))
    if xscale is None:
        if ax is None:
            xscale = "linear"
        else:
            xscale = ax.get_xscale()
    if yscale is None:
        if ax is None:
            yscale = "linear"
        else:
            yscale = ax.get_yscale()
    if ax is None:
        ax = fig.add_subplot(111)
    if plot_used:
        used = y > 0
        xmid = xmid[used]
        y = y[used]
        xerr = xerr[used]
        yerr = yerr[used]
    if noerr:
        ax.plot(xmid, y, lw=lw, label=label, **kwargs)
    else:
        ax.errorbar(
            xmid, y, yerr=yerr, xerr=xerr, fmt=".", lw=lw, label=label, **kwargs
        )
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    if plot_counts:
        ylabel = "Counts (counts/{0})"
    else:
        ylabel = "Count Rate (counts/s/{0})"
    ax.set_ylabel(ylabel.format(yunit), fontsize=fontsize)
    ax.tick_params(reset=True, labelsize=fontsize)
    return fig, ax, ebins
