def grid_overlay(axes, grid_spacing):
    """
    Create a heliographic overlay using wcsaxes.
    Also draw a grid and label the top axes.

    Parameters
    ----------
    axes : `~astropy.visualization.wcsaxes.WCSAxes` object.
        The `~astropy.visualization.wcsaxes.WCSAxes` object to create the HGS overlay on.

    grid_spacing: `~astropy.units.Quantity`
        Spacing for longitude and latitude grid in degrees.

    Returns
    -------
    overlay : `~astropy.visualization.wcsaxes.WCSAxes` overlay
        The overlay object.

    """
    lon_space = lat_space = grid_spacing
    overlay = axes.get_coords_overlay('heliographic_stonyhurst')
    lon = overlay[0]
    lat = overlay[1]
    lon.coord_wrap = 180
    lon.set_major_formatter('dd')
    lon.set_ticks_position('tr')
    lat.set_ticks_position('tr')
    grid_kw = {'color': 'white', 'zorder': 100, 'alpha': 0.5}#, linestyle: 'dashed', linewidth: 0.1}
    lon.set_ticks(spacing=lon_space, color=grid_kw['color'])
    lat.set_ticks(spacing=lat_space, color=grid_kw['color'])
    overlay.grid(**grid_kw, linestyle='dashed', linewidth=0.1)
    return overlay