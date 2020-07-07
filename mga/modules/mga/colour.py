#!/usr/bin/env python

""" MultiQC Multi Genome Alignment module """

class Colour(object):
    
    @staticmethod
    def fromBytes(r, g, b):
        return Colour(r, g, b)
    
    def __init__(self, r, g, b):
        self.red = r
        self.green = g
        self.blue = b
        
    def applyAlpha(self, alpha):
        return Colour(self._alphaValue(self.red, alpha),
                      self._alphaValue(self.green, alpha),
                      self._alphaValue(self.blue, alpha))
    
    def toHtml(self):
        return "#{:02x}{:02x}{:02x}".format(self.red, self.green, self.blue)
    
    def _alphaValue(self, c, a):
        # 255 because we're always considering a white background.
        bgc = 255
        return int(a * c + (1.0 - a) * bgc)
