#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2015, ENPC
#     Author(s): Sylvain Doré
#
# This file is part of the air quality modeling system Polyphemus.
#
# Polyphemus is developed in the INRIA project-team CLIME and in
# the ENPC - EDF R&D joint laboratory CEREA.
#
# Polyphemus is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Polyphemus is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# more details.
#
# For more information, visit the Polyphemus web site:
#      http://cerea.enpc.fr/polyphemus/

import os.path
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from numpy import pi, cos, exp


def show_img(path, **kwargs):
    try:
        is_notebook = get_ipython().kernel is not None
    except:
        is_notebook = False

    if is_notebook:
        from IPython.display import Image, display
        display(Image(path))
    else:
        import webbrowser
        webbrowser.open("file://" + os.path.abspath(path))


def contourf_map(ax, basemap, data, **kwargs):
    x_range, y_range = basemap.makegrid(data.shape[1], data.shape[0])
    cs = basemap.contourf(x_range, y_range, data, ax=ax, **kwargs)
    basemap.drawcountries(ax=ax)
    basemap.drawcoastlines(ax=ax)
    return cs


def plot_array_list(array_list,
                    nb_image_t=1,
                    nb_image_z=1,
                    title=None,
                    array_name_list=None,
                    plot_func=None,
                    png_filename=None,
                    show_figure=True):
    if plot_func == None:
        def default_plot_func(fig, ax, img, n, t, z):
            cs = ax.contourf(img)
            fig.colorbar(cs, ax=ax)

        plot_func = default_plot_func

    ncols = len(array_list)
    nrows = nb_image_t * nb_image_z
    fig, ax = plt.subplots(nrows, ncols)

    if array_name_list:
        array_name_title = ','.join(array_name_list)
        if title is None:
            title = array_name_title
        else:
            title += "\n\n(" + array_name_title + ")"

    for c in xrange(0, ncols):
        data = array_list[c]
        Nt = data.shape[0]
        Nz = data.shape[1]
        sampling_range_t = np.linspace(0, Nt - 1, nb_image_t, dtype=np.int)
        sampling_range_z = np.linspace(0, Nz - 1, nb_image_z, dtype=np.int)
        img = [dict(img=data[t][z],
                    n=c,
                    t=t,
                    z=z) for z in sampling_range_z for t in sampling_range_t]

        label = ""
        if Nt > 1 and Nz > 1:
            label = "t={t}, z={z}"
        elif Nt > 1:
            label = "t={t}"
        elif Nz > 1:
            label = "z={z}"

        for r in xrange(0, nrows):
            ax_rc=ax[r][c]
            image = img[r]
            plt.sca(ax_rc)
            plot_func(fig, ax_rc, **image)
            if label:
                ax_rc.set_title(label.format(**image), fontsize='xx-large')

    image_shape = array_list[0].shape[2:4]
    image_ratio = float(image_shape[0]) / image_shape[1]
    image_width = 20.
    fig.set_size_inches(image_width,
                        image_width * image_ratio * (float(nrows) / ncols))

    plt.tight_layout()

    if title:
        title_height=2.5
        w, h = fig.get_size_inches()
        h += title_height
        fig.set_size_inches(w, h)
        plt.subplots_adjust(top=(h - title_height)/h)
        fig.suptitle(title, fontsize=25, y=(h - title_height/3.)/h)


    if png_filename:
        plt.savefig(png_filename, dpi=64)

    if show_figure:
        # Must come in last position since it resets the current figure.
        plt.show()

    plt.close()


def normalize_colormap(data_array_list,
                       cmap=mpl.cm.get_cmap("seismic"),
                       zero_centered=True,
                       log_normed=False):
    # The computation time is optimized for a sequence of numpy arrays
    # or a single numpy array.

    # To manage the case where a single numpy array was given:
    if type(data_array_list) is np.ndarray:
        data_array_list = [data_array_list]

    # Could have an auto log normalisation by comparing the absolute extrema
    # with the absolute mean.

    vmin = min([np.min(a) for a in data_array_list])
    vmax = max([np.max(a) for a in data_array_list])

    if zero_centered:
        vmax = max(abs(vmin), vmax)
        vmin = -vmax

    norm = None
    if log_normed:
        if vmin is None:
            norm = mpl.colors.LogNorm()
        else:
            if vmin < 0:
                vmean = np.mean([np.mean(np.abs(a)) for a in data_array_list])
            else:
                vmean = np.mean([np.mean(a) for a in data_array_list])

            # If the data crosses or comes close to zero, then a linearization
            # of the logarithm is needed for the smallest values.
            if vmin < 0 or vmin < vmean / 100.:
                norm = mpl.colors.SymLogNorm(vmin=vmin,
                                             vmax=vmax,
                                             linthresh=vmean / 100.)
            else:
                norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax)

    plot_kwargs = dict(vmin=vmin, vmax=vmax, norm=norm, cmap=cmap)
    return plot_kwargs


def dummy_data(Nx=64,
               Ny=64,
               Nt=32,
               Nz=1,
               amplitude=1,
               speed=0.1,
               dtype=np.float32):
    L = 5  # Size of the zero centered square domain.
    x = np.linspace(-L / 2., L / 2., Nx)
    y = np.linspace(-L / 2., L / 2., Ny)
    y = y[:, np.newaxis]  # transpose

    # A dummy function to generate data.
    def f(z, t):
        A = amplitude / z ** 2
        return A * cos(pi * x * y * t) ** 2 * exp(-(x ** 2 + y ** 2) / 2.)

    dz = 10. / Nz
    return np.array([[f(z, t) for z in np.arange(1., Nz * dz + 1., dz)]
                     for t in np.arange(0., Nt * speed * dz, speed * dz)],
                    dtype)
