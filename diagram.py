# coding=utf-8

from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram._LinearDrawer import LinearDrawer
from Bio.Graphics.GenomeDiagram._CircularDrawer import CircularDrawer as _Cir
from Bio.Graphics.GenomeDiagram._CircularDrawer import _stroke_and_fill_colors, ArcPath, pi, cos, sin, draw_polygon
from reportlab.graphics.charts.legends import Legend


class CircularDrawer(_Cir):
    def draw_legend(self, name_pairs):
        legend = Legend()
        legend.colorNamePairs = name_pairs
        fontName = 'Helvetica'
        fontSize = 7
        legend.fontName = fontName
        legend.fontSize = fontSize
        legend.alignment = 'right'
        legend.dxTextSpace = 7
        legend.boxAnchor = 'nw'
        legend.subCols.dx = 0
        legend.subCols.dy = -2
        legend.subCols.rpad = 0
        legend.columnMaximum = len(name_pairs)
        legend.deltax = 1
        legend.deltay = 0
        legend.dy = 5
        legend.dx = 5

        radius = self.track_radii[self.current_track_level][1]
        legend.x = self.xcenter - radius / 4
        legend.y = self.ycenter + radius / 4 * 3
        self.drawing.add(legend)

    def draw_feature(self, feature):
        """Return list of feature elements and list of labels for them."""
        feature_elements = []  # Holds drawable elements for a single feature
        label_elements = []  # Holds labels for a single feature

        if feature.hide:  # Don't show feature: return early
            return feature_elements, label_elements

        start, end = self._current_track_start_end()
        # A single feature may be split into subfeatures, so loop over them
        for locstart, locend in feature.locations:
            if locend < start:
                continue
            locstart = max(locstart, start)
            if end < locstart:
                continue
            locend = min(locend, end)
            # Get sigil for the feature/ each subfeature
            kwargs = feature.__dict__.copy()
            [kwargs.__delitem__(k) for k in ('color', 'border', 'strand')]
            feature_sigil, label = self.get_feature_sigil(feature, locstart, locend, **kwargs)
            feature_elements.append(feature_sigil)
            if label is not None:  # If there's a label
                label_elements.append(label)

        return feature_elements, label_elements

    def _calc_center(self, inner_radius, startangle, endangle, tilted_angle):
        mid_angle = (startangle + endangle) / 2
        mid_sin, mid_cos = sin(mid_angle), cos(mid_angle)
        xc, yc = self.xcenter, self.ycenter
        txc, tyc = (xc + inner_radius * mid_sin, yc + inner_radius * mid_cos)
        tilt_sin, tilt_cos = sin(mid_angle + tilted_angle), cos(mid_angle + tilted_angle)
        tx0, ty0 = (txc - inner_radius * tilt_sin, tyc - inner_radius * tilt_cos)
        return tx0, ty0, startangle + tilted_angle, endangle + tilted_angle

    def _draw_arc(self, inner_radius, outer_radius, startangle, endangle,
                  color, border=None, colour=None, **kwargs):
        """Return closed path describing an arc box.

        Arguments:
         - inner_radius  Float distance of inside of arc from drawing center
         - outer_radius  Float distance of outside of arc from drawing center
         - startangle    Float angle subtended by start of arc at drawing center
           (in radians)
         - endangle      Float angle subtended by end of arc at drawing center
           (in radians)
         - color        colors.Color object for arc (overridden by backwards
           compatible argument with UK spelling, colour).

        Returns a closed path object describing an arced box corresponding to
        the passed values.  For very small angles, a simple four sided
        polygon is used.
        """
        # Let the UK spelling (colour) override the USA spelling (color)
        if colour is not None:
            color = colour

        strokecolor, color = _stroke_and_fill_colors(color, border)

        if abs(float(endangle - startangle)) > .01:
            # Wide arc, must use full curves
            p = ArcPath(strokeColor=strokecolor,
                        fillColor=color,
                        strokewidth=0)
            # Note reportlab counts angles anti-clockwise from the horizontal
            # (as in mathematics, e.g. complex numbers and polar coordinates)
            # but we use clockwise from the vertical.  Also reportlab uses
            # degrees, but we use radians.
            p.addArc(self.xcenter, self.ycenter, inner_radius,
                     90 - (endangle * 180 / pi), 90 - (startangle * 180 / pi),
                     moveTo=True)
            p.addArc(self.xcenter, self.ycenter, outer_radius,
                     90 - (endangle * 180 / pi), 90 - (startangle * 180 / pi),
                     reverse=True)
            p.closePath()
            return p
        else:
            # Cheat and just use a four sided polygon.

            x0, y0 = self.xcenter, self.ycenter  # origin of the circle
            if kwargs.has_key('tilted_angle'):
                tilted_angle = kwargs['tilted_angle']
                tilted_angle = tilted_angle * pi / 180
                x0, y0, startangle, endangle = self._calc_center(inner_radius, startangle, endangle, tilted_angle)

            # Calculate trig values for angle and coordinates
            startcos, startsin = cos(startangle), sin(startangle)
            endcos, endsin = cos(endangle), sin(endangle)

            x1, y1 = (x0 + inner_radius * startsin, y0 + inner_radius * startcos)
            x2, y2 = (x0 + inner_radius * endsin, y0 + inner_radius * endcos)
            x3, y3 = (x0 + outer_radius * endsin, y0 + outer_radius * endcos)
            x4, y4 = (x0 + outer_radius * startsin, y0 + outer_radius * startcos)
            return draw_polygon([(x1, y1), (x2, y2), (x3, y3), (x4, y4)], color, border)


class Diagram(GenomeDiagram.Diagram):
    def draw(self, format=None, pagesize=None, orientation=None,
             x=None, y=None, xl=None, xr=None, yt=None, yb=None,
             start=None, end=None, tracklines=None, fragments=None,
             fragment_size=None, track_size=None, circular=None,
             circle_core=None, cross_track_links=None, color_name_pairs=None):
        """Draw the diagram, with passed parameters overriding existing attributes.

        gdd.draw(format='circular')
        """
        # Pass the parameters to the drawer objects that will build the
        # diagrams.  At the moment, we detect overrides with an or in the
        # Instantiation arguments, but I suspect there's a neater way to do
        # this.
        if format == 'linear':
            drawer = LinearDrawer(self, pagesize or self.pagesize,
                                  orientation or self.orientation,
                                  x or self.x, y or self.y, xl or self.xl,
                                  xr or self.xr, yt or self.yt,
                                  yb or self.yb, start or self.start,
                                  end or self.end,
                                  tracklines or self.tracklines,
                                  fragments or self.fragments,
                                  fragment_size or self.fragment_size,
                                  track_size or self.track_size,
                                  cross_track_links or self.cross_track_links)
        else:
            drawer = CircularDrawer(self, pagesize or self.pagesize,
                                    orientation or self.orientation,
                                    x or self.x, y or self.y, xl or self.xl,
                                    xr or self.xr, yt or self.yt,
                                    yb or self.yb, start or self.start,
                                    end or self.end,
                                    tracklines or self.tracklines,
                                    track_size or self.track_size,
                                    circular or self.circular,
                                    circle_core or self.circle_core,
                                    cross_track_links or self.cross_track_links)
        drawer.draw()  # Tell the drawer to complete the drawing
        if color_name_pairs is not None:
            drawer.draw_legend(color_name_pairs)
        self.drawing = drawer.drawing  # Get the completed drawing
