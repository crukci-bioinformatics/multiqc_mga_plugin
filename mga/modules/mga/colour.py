#!/usr/bin/env python

""" MultiQC Multi Genome Alignment module """

import colorsys

class Colour(object):
    
    @staticmethod
    def fromBytes(r, g, b):
        return Colour(r / 255.0, g / 255.0, b / 255.0)
    
    def __init__(self, r, g, b):
        self.red = r
        self.green = g
        self.blue = b
        
    def applyAlpha(self, alpha):
        h, l, s = colorsys.rgb_to_hls(self.red, self.green, self.blue)
        rgb = colorsys.hls_to_rgb(h, l, s * alpha)
        return Colour(rgb[0], rgb[1], rgb[2])
    
    def toHtml(self):
        return "#{:02x}{:02x}{:02x}".format(int(self.red * 255), int(self.green * 255), int(self.blue * 255))
