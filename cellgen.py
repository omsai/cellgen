"""Generate test cell images with combinations of centromeres.

Each centromere may or may not have an overlapping ectopic centromere.
The amount of overlap varies by image.  The nucleus size is uniformly
set such that it encloses the largest collection of centromeres."""

import itertools
import os
import subprocess
import sys
import unittest

import PIL
import PIL.ImageDraw
import numpy as np
import scipy.optimize


class PointsUniform(object):
    """Generates uniformly distributed array of points in 2 dimensions.

    All points in the array are between 0 and 1.

    >>> pu = PointsUniform(100)
    >>> pu
    <PointsUniform object with 100 points at least 0.0888388838748 apart>
    >>> pu.dist
    0.0888388838748
    >>> pu.points[:5]
    [[ 0.815431  0.064958]
     [ 0.073785  0.54831 ]
     [ 0.800585  0.641348]
     [ 0.20884   0.111673]
     [ 0.394668  0.29241 ]]
    """
    # pylint: disable=too-few-public-methods

    # Path to the `fast_delaunay` executable.
    _path_fast_delunay = "./fast-poisson-disk/fast_delaunay"

    def __init__(self, n_points, debug=False):
        """Maximize distances using Golden method of optimization."""
        # maximum `dist` between points should be in range (0, 1]
        self.n_points = n_points
        # Optimize for largest distance.
        kwargs = {}
        if debug:
            kwargs["disp"] = 3
        self.dist = scipy.optimize.fminbound(self.generate_points, 0, 1,
                                             **kwargs)
        # Regenerate points with the optimum value.
        self.generate_points(self.dist)

    def generate_points(self, dist):
        """Set numpy array from executable program.

        Return distance optimization score."""
        cmd = "{prog} {n} {radius}".format(prog=self._path_fast_delunay,
                                           n=self.n_points,
                                           radius=dist)
        # Use devnull to sink "Done." output.  Better would be to patch
        # the program not to write "Done." to stderr.
        with open(os.devnull, 'w') as devnull:
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                       stderr=devnull, shell=True)
        ret = process.communicate()[0]
        self.points = np.fromstring(ret, sep=' ').reshape(-1, 2)
        return -dist + abs(self.points.shape[0] - self.n_points)

    def __repr__(self):
        return ("<{class_} object with {n} points "
                "at least {dist} apart>").format(
                    class_=self.__class__.__name__,
                    n=self.n_points,
                    dist=self.dist)


class CellImages(object):
    """Generate cell nuclei with centromeres."""
    # pylint: disable=too-few-public-methods
    # pylint: disable=too-many-instance-attributes

    def __init__(self, n_cen, n_ect, im_cen_overlap, width_cen_px=5):
        """Set nuclei, centromeres and overlapping ectopic centromeres.

        Parameters
        ----------
        n_cen : integer
            Number of centromeres in cell.
        n_ect : integer
            Number of ectopic centromeres in cell.
        im_cen_overlap : array-like
            Amount of overlap of ectopic centromeres with regular
            centromeres in range [-1, 1].
        width_cen_px : integer, optional
            Pixel width of centromeres in generated images.

        """
        self.nuclei_children = np.array(
            list(
                itertools.product(range(1, 1 + n_cen),
                                  range(1, 1 + n_ect))))
        if len(self.nuclei_children) < 1:
            raise ValueError, ("Not enough centromeres to produce nuclei."
                               "  Try increasing the value of `n_cen` or `n_ect`.")
        # Calculate unitless point distributions.
        self._pu_nuc = PointsUniform(len(self.nuclei_children))
        self._pu_cen = PointsUniform(max(n_cen, n_ect))

        # Calculate geometric distances in pixels.
        self._padding_nuc_px = int(
            np.ceil(
                2 * width_cen_px
                + (1 - np.min(im_cen_overlap))
                * width_cen_px))
        self.width_nuc_px = (
            int(
                np.ceil(
                    2 * width_cen_px
                    / self._pu_cen.dist))
            + self._padding_nuc_px)
        self._padding_im_px = 2 * self.width_nuc_px
        self.width_im_px = (
            int(
                np.ceil(
                    2 * self.width_nuc_px
                    / self._pu_nuc.dist))
            + self._padding_im_px)

        self.n_cen = n_cen
        self.n_ect = n_ect
        self.im_cen_overlap = im_cen_overlap
        self.width_cen_px = width_cen_px

    @property
    def frames(self):
        """Total number of image frames"""
        return len(self.im_cen_overlap)

    def image(self, frame, channel):
        """Get single image

        Parameters
        ----------
        frame : integer
            Frame of cells with constant amount of overlap.
        channel : integer
            0 = nuc, 1 = cen, 2 = ect.

        """
        # pylint: disable=too-many-locals
        im = PIL.Image.new('L',   # pylint: disable=invalid-name
                           (self.width_im_px, self.width_im_px))
        draw = PIL.ImageDraw.Draw(im)
        for i, nuc in enumerate(self._pu_nuc.points
                                * (self.width_im_px - self._padding_im_px)
                                + self._padding_im_px / 2):
            if channel == 0:
                draw.ellipse(bounding_box(nuc, self.width_nuc_px
                                          + self._padding_nuc_px), fill=128)
            else:
                bbox = lambda x: bounding_box(x, self.width_cen_px)
                n_cen = self.nuclei_children[i, 0] # pylint: disable=invalid-sequence-index
                n_ect = self.nuclei_children[i, 1] # pylint: disable=invalid-sequence-index
                for j, cen in enumerate(self._pu_cen.points
                                        * (self.width_nuc_px
                                           - self._padding_nuc_px)
                                        - self._padding_nuc_px / 2
                                        + nuc):
                    if channel == 1 and j < n_cen:
                        draw.ellipse(bbox(cen), fill=128)
                    elif channel == 2 and j < n_ect:
                        overlap = self.im_cen_overlap[frame]
                        dist = self.width_cen_px * (1 - overlap)
                        ect_diff = point_along_line([0, cen[0] - nuc[0]],
                                                    [0, cen[1] - nuc[1]], dist)
                        ect = ect_diff + np.array(nuc)
                        draw.ellipse(bbox(ect), fill=128)
        return im

    def save(self, path="./",
             format_="{frame}-ch{channel}_{channel_name}-ov{overlap}.tif"):
        """Write all images."""
        channel_names = ["nuc", "cen", "ect"]
        for frame in range(self.frames):
            for channel, channel_name in enumerate(channel_names):
                im_name = format_.format(frame=frame, channel=channel,
                                         channel_name=channel_name,
                                         overlap=self.im_cen_overlap[frame])
                self.image(frame, channel).save(os.path.join(path, im_name))

    def __repr__(self):
        return ("<{class_} {f} image frames "
                "x {n_nuc} nuclei "
                "x (up to) {n_cen} centromeres, "
                "(up to) {n_ect} overlapping ectopic centromeres>").format(
                    class_=self.__class__.__name__,
                    f=self.frames,
                    n_nuc=len(self.nuclei_children),
                    n_cen=self.n_cen,
                    n_ect=self.n_ect)


class TestCell(unittest.TestCase):
    """Test the Cell class"""
    # pylint: disable=missing-docstring

    def test_error_too_few_cen(self):
        with self.assertRaises(ValueError):
            CellImages(n_cen=0, n_ect=0, im_cen_overlap=[-1, 0, 1])

    def test_px_distances_are_sane(self):
        # pylint: disable=protected-access
        cells = CellImages(2, 2, im_cen_overlap=[-1, 0, 1])
        self.assertLess(cells.width_nuc_px, cells.width_im_px)
        self.assertLess(cells._padding_nuc_px, cells.width_nuc_px)
        self.assertLess(cells._padding_im_px, cells.width_im_px)

    # def test_enumerated_centromeres(self):
    #     raise NotImplemented

    # def test_enumerated_cells(self):
    #     raise NotImplemented

    # def test_overlap_none(self):
    #     raise NotImplemented


def bounding_box(point, width):
    """Return bounding box of upper-left and lower-right corners."""
    # pylint: disable=invalid-name
    d = width / 2.
    x, y = point
    return x - d, y - d, x + d, y + d

def point_along_line(x, y, dist):
    """Get point at specified distance along line.

    Returns [x, y] coordinate along, but not in, the line.

    """
    # pylint: disable=invalid-name
    # pylint: disable=no-member
    slope, _ = np.polyfit(x, y, 1)
    theta = np.arctan(slope)
    x_ = x[-1] + dist * np.sin(theta) * np.sign(y[-1])
    y_ = y[-1] + dist * np.cos(theta) * np.sign(y[-1])
    return np.array([x_, y_])

if __name__ == '__main__':
    PATH = sys.argv[1:]
    if len(PATH) != 1:
        PATH = "./"
    else:
        PATH = PATH[0]
    CellImages(4, 4, np.linspace(-1, 1, 5)).save(path=PATH)
