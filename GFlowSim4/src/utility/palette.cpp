#include "palette.hpp"
// Other files
#include "vectormath.hpp"
#include <algorithm> // For swap

namespace GFlowSimulation {

  Palette::Palette(int width, int height) {
    if (width<=0 || height<=0) throw BadPalette();
    // Set the bounds
    owned[0] = bounds[0] = 0;
    owned[1] = bounds[1] = width;
    owned[2] = bounds[2] = 0;
    owned[3] = bounds[3] = height;
    // Create the image
    image = new BMP;
    image->SetSize(width, height);
    pixelData = new Array<PixL, 2>(width, height);
    // References
    refs = new int(1);
    // Aspect ratio
    aspect_ratio = 1;
    combinationRule = 0; // <-
  }

  Palette::Palette(const Palette& p) {
    if (p.getWidth()<=0 || p.getHeight()<=0) throw BadPalette();
    // Get bounds
    for (int i=0; i<4; ++i) {
      owned[i] = p.owned[i];
      bounds[i] = p.bounds[i];
    }
    // Image and refs
    image = p.image;
    pixelData = p.pixelData; // <-
    refs = p.refs;
    ++refs[0];
    // Aspect ratio
    aspect_ratio = p.aspect_ratio;
    combinationRule = p.combinationRule;
  }

  Palette::Palette(const Palette&& p) {
    if (p.getWidth()<=0 || p.getHeight()<=0) throw BadPalette();
    // Get bounds
    for (int i=0; i<4; ++i) {
      owned[i] = p.owned[i];
      bounds[i] = p.bounds[i];
    }
    // Image and refs
    image = p.image;
    pixelData = p.pixelData; // <-
    refs = p.refs;
    ++refs[0];
    // Aspect ratio
    aspect_ratio = p.aspect_ratio;
    combinationRule = p.combinationRule;
  }
    
  Palette::~Palette() {
    if (!refs) return;
    // Clean up if we are the last palette to reference the object
    if (refs[0]==1) {
      --refs[0];
      if (image) delete image;
      if (pixelData) delete pixelData; // <-
      delete refs;
    }
    else --refs[0];
  }

  Palette& Palette::operator=(const Palette&& p) {
    if (p.getWidth()<=0 || p.getHeight()<=0) throw BadPalette();
    // Get bounds
    for (int i=0; i<4; ++i) {
      owned[i] = p.owned[i];
      bounds[i] = p.bounds[i];
    }
    // Image and refs
    image = p.image;
    pixelData = p.pixelData; // <-
    combinationRule = p.combinationRule;
    refs = p.refs;
    ++refs[0];
    // Aspect ratio
    aspect_ratio = p.aspect_ratio;
    // Return
    return *this;
  }

  void Palette::writeToFile(string fileName) {
    pixelsToImage(image, 0, bounds[1], 0, bounds[3]);
    image->WriteToFile(fileName.c_str());
  }

  Palette Palette::getSubPalette(int mx, int Mx, int my, int My) {
    // Check for errors
    if (Mx<mx || My<my) throw BadPalette();
    // Correct
    if (mx<0) mx=0;
    else if (mx>=getWidth()) mx = getWidth();
    if (my<0) my = 0;
    else if (my>getHeight()) my = getHeight();
    // Return a palatte
    return Palette(owned[0]+mx, owned[0]+Mx, owned[2]+my, owned[2]+My, image, refs, bounds, aspect_ratio, pixelData);
  }

  Palette Palette::getSubPalette(float mfx, float Mfx, float mfy, float Mfy) {
    return getSubPalette(
      static_cast<int>(mfx*getWidth()), 
      static_cast<int>(Mfx*getWidth()), 
      static_cast<int>(mfy*getHeight()), 
      static_cast<int>(Mfy*getHeight())
    );
  }

  Palette Palette::getSubPaletteCentered(float fx) {
    if (fx>=0.5) throw BadPalette();
    int pix = getWidth()*fx;
    return getSubPalette(pix, getWidth()-pix, pix, getHeight()-pix);
  }

  Palette Palette::getMaxCenteredSubPalette(float rx, float ry) {
    if (rx<=0 || ry<= 0) throw BadPalette();
    double native_ratio = getWidth()/getHeight();
    double requested_ratio = rx/ry;

    if (native_ratio<requested_ratio) { // Fit to x. There will be extra y-space. 
      int width = getWidth();
      int height = width/requested_ratio;
      int gap = (getHeight() - height)/2;
      return getSubPalette(0, getWidth(), gap, getHeight()-gap);
    }
    else { // Fit to y. There will be extra x-space.
      int height = getHeight();
      int width = requested_ratio*height;
      int gap = (getWidth() - width)/2;
      return getSubPalette(gap, getWidth()-gap, 0, getHeight());
    }
  }

  void Palette::drawCircleByFactors(float xf, float yf, float rfx, ColorFunction colorF, bool wrap) {
    // Width of owned region
    int wx = getWidth();
    int wy = getHeight();
    float iwx = 1.f/wx, iwy = 1.f/wy;
    // Find the center (in pixel coordinates)
    int px = xf*wx;
    int py = yf*wy;
    // Find dp's (radius in pixel coordinates)
    int dpx = ceil(rfx*wx);
    int dpy = aspect_ratio*dpx;

    // Radius factor squared
    float rsqr = sqr(rfx);
    bool doColor = true;
    // Go through pixels
    for (int dy=-dpy; dy<=dpy; ++dy) {
      for (int dx=-dpx; dx<=dpx; ++dx) {
        if (sqr(iwx*dx) + sqr(iwy*dy) <= rsqr) {
          // The potential pixel
          int w_px = px + dx, w_py = py + dy;
          // Get the pixel factor
          std::pair<float, float> pF = pixelFactor(w_px, w_py);
          float pfx = (pF.first - xf)/rfx;
          float pfy = (pF.second - yf)/rfx;
          // Wrap pixel
          if (wrap) {
            w_px = w_px >= wx ? w_px-wx : w_px;
            w_px = w_px <  0 ? w_px+wx : w_px;
            w_py = w_py >= wy ? w_py-wy : w_py;
            w_py = w_py <  0 ? w_py+wy : w_py;
            // Set pixel
            RGBApixel color = colorF(pfx, pfy, doColor);
            if (doColor) setPixel(w_px, w_py, color);
          }
          else if (0<=w_px && w_px<wx && 0<=w_py && w_py<wy) {
            // Set pixel
            RGBApixel color = colorF(pfx, pfy, doColor);
            if (doColor) setPixel(w_px, w_py, color);
          }
        }
      }
    }
  }

  void Palette::drawLineByFactors(float xf0, float yf0, float xf1, float yf1, ColorFunction colorF, float w) {
    // The algorithm takes coordinates where pixels have integer coordinates. Convert to that basis.
    int x0 = xf0*getWidth(),  x1 = xf1*getWidth();
    int y0 = yf0*getHeight(), y1 = yf1*getHeight();
    // Draw a line using Wu's algorithm
    if (colorF!=nullptr) {
      bool bl;
      drawLine_WuAlgorithm(x0, y0, x1, y1, colorF(0,0,bl));
    }
    else drawLine_WuAlgorithm(x0, y0, x1, y1, RGB_Blue);
  }

  void Palette::coverPalette(const RGBApixel color) {
    for (int y=owned[2]; y<owned[3]; ++y)
      for (int x=owned[0]; x<owned[1]; ++x) 
        image->SetPixel(x, y, color);
  }

  void Palette::drawGraph2d(vector<pair<float, float> >& data, GraphOptions& options) {
    if (options.paintBackground) coverPalette(RGB_Black);
    Palette face = getSubPalette(0.01f, 0.99f, 0.02f, 0.98f);
    face.drawGraphData2d(data, options);
  }

  int Palette::getWidth() const {
    return owned[1] - owned[0];
  }

  int Palette::getHeight() const {
    return owned[3] - owned[2];
  }

  Palette::Palette(int mx, int Mx, int my, int My, BMP *img, int *rfs, int *bnds, double aratio, Array<PixL, 2> *pxlData) {
    if (Mx<=mx || My<=my) throw BadPalette();
    // Set sub-bounds
    owned[0] = mx;
    owned[1] = Mx;
    owned[2] = my;
    owned[3] = My;
    // Copy full bounds
    for (int i=0; i<4; ++i) bounds[i] = bnds[i];
    // Image and refs
    image = img;
    pixelData = pxlData; // <-
    combinationRule = 0;
    refs = rfs;
    ++refs[0];
    // Aspect ratio
    aspect_ratio = aratio;
  }

  inline std::pair<float, float> Palette::pixelFactor(int x, int y) const {
    // Create a pair to set and return
    std::pair<float, float> pf;
    // We need to find where the centers of pixels are
    int X = owned[1]-owned[0]+1;
    int Y = owned[3]-owned[2]+1;
    // Set pixel factor
    pf.first  = static_cast<float>(x)/X;
    pf.second = static_cast<float>(y)/Y;
    // Return
    return pf;
  }

  inline float Palette::pixelFactorX(int x) const {
    int X = owned[1]-owned[0]+1;
    return static_cast<float>(x)/X;
  }

  inline float Palette::pixelFactorY(int y) const {
    int Y = owned[3]-owned[2]+1;
    return static_cast<float>(y)/Y;
  }

  //! @brief Set a pixel relative to the owned coordinates (origin is {owned[0], owned[2]})
  //!
  //! EasyBMP has the top left corner as its origin, it would be easier to use the bottom left, so
  //! this function also takes care of that.
  inline void Palette::setPixel(int x, int y, RGBApixel color, float depth) {
    int X = x+owned[0], Y = owned[3]-y-1;
    // Check the bounds
    if (X<0 || bounds[1]<=X || Y<0 || bounds[3]<=Y) return;
    // Set the pixel
    else {
      // image->SetPixel(X, Y, color);
      pixelData->at(X, Y).set(color, depth);
    }
  }

  inline void drawLine_BresenhamAlgorithm(int x0, int y0, int x1, int y1, RGBApixel color) {
    throw Unimplemented();
  }

  inline void Palette::drawLine_WuAlgorithm(int x0, int y0, int x1, int y1, RGBApixel color) {
    // Is the line steep?
    bool isSteep = fabs(y1-y0) > fabs(x1-x0);

    if (isSteep) {
      std::swap(x0, y0);
      std::swap(x1, y1);
    }
    if (x0>x1) {
      std::swap(x0, x1);
      std::swap(y0, y1);
    }

    float dx = x1 - x0;
    float dy = y1 - y0;
    float gradient = dx==0 ? 1.0 : dy/dx;

    // Some helpful lambdas
    auto fpart = [] (float x) -> float {
      return x - floor(x);
    };
    auto rfpart = [] (float x) -> float {
      return 1. - x + floor(x);
    };

    // handle first endpoint
    float xend = round(x0);
    float yend = y0 + gradient * (xend - x0);
    float xgap = rfpart(x0 + 0.5);
    float xpxl1 = xend; // this will be used in the main loop
    float ypxl1 = floor(yend);
    if (isSteep) {
      setPixel(ypxl1,   xpxl1, rfpart(yend) * xgap * color);
      setPixel(ypxl1+1, xpxl1,  fpart(yend) * xgap * color);
    }
    else {
      setPixel(xpxl1, ypxl1  , rfpart(yend) * xgap * color);
      setPixel(xpxl1, ypxl1+1,  fpart(yend) * xgap * color);
    }

    float intery = yend + gradient; // first y-intersection for the main loop
    
    // handle second endpoint
    xend = round(x1);
    yend = y1 + gradient * (xend - x1);
    xgap = fpart(x1 + 0.5);
    float xpxl2 = xend; //this will be used in the main loop
    float ypxl2 = floor(yend);
    if (isSteep) {
      setPixel(ypxl2  , xpxl2, rfpart(yend) * xgap * color);
      setPixel(ypxl2+1, xpxl2,  fpart(yend) * xgap * color);
    }
    else {
      setPixel(xpxl2, ypxl2,  rfpart(yend) * xgap * color);
      setPixel(xpxl2, ypxl2+1, fpart(yend) * xgap * color);
    }
    
    // main loop
    if (isSteep) {
      for (int x=xpxl1 + 1; x<=xpxl2 - 1; ++x) {
        setPixel(floor(intery)  , x, rfpart(intery) * color);
        setPixel(floor(intery)+1, x,  fpart(intery) * color);
        intery += gradient;
      }
    }
    else {
      for (int x=xpxl1 + 1; x<=xpxl2 - 1; ++x) {
        setPixel(x, floor(intery),  rfpart(intery) * color);
        setPixel(x, floor(intery)+1, fpart(intery) * color);
        intery += gradient;
      }
    }
    // Thanks, Wikipedia!
  }

  inline void Palette::drawCircle_MidpointAlgorithm(int x0, int y0, int radius, RGBApixel color) {
    int x = radius-1;
    int y = 0;
    int dx = 1;
    int dy = 1;
    int err = dx - (radius << 1);

    while (x >= y) {
      setPixel(x0 + x, y0 + y, color);
      setPixel(x0 + y, y0 + x, color);
      setPixel(x0 - y, y0 + x, color);
      setPixel(x0 - x, y0 + y, color);
      setPixel(x0 - x, y0 - y, color);
      setPixel(x0 - y, y0 - x, color);
      setPixel(x0 + y, y0 - x, color);
      setPixel(x0 + x, y0 - y, color);

      if (err <= 0) {
        y++;
        err += dy;
        dy += 2;
      }
      
      if (err > 0) {
        x--;
        dx += 2;
        err += dx - (radius << 1);
      }
    }
    // Thanks, Wikipedia! <https://en.wikipedia.org/wiki/Midpoint_circle_algorithm>
  }

  inline void Palette::drawGraphData2d(vector<pair<float, float> >& data, GraphOptions& options) {
    // Check if there is anything to draw and if it is possible to draw.
    if (image==nullptr || data.empty()) return;
    // Find bounds
    float mx, Mx, my, My;
    pair<float, float> first = *data.begin();

    // Set options
    RGBApixel background = RGB_White, lineColor = RGB_Black;
    if (options.useLineColor) lineColor = options.lineColor;
    // Check bounds
    if (!options.hasBounds) { 
      // Find the min/max bounds
      mx = Mx = first.first;
      my = My = first.second;
      for (auto p : data) {
        if (p.first>Mx)  Mx = p.first;
        if (p.first<mx)  mx = p.first;
        if (p.second>My) My = p.second;
        if (p.second<my) my = p.second;
      }

      // If min/max specified
      if (options.useMinX) mx = options.minX;
      if (options.useMaxX) Mx = options.maxX;
      if (options.useMinY) my = options.minY;
      if (options.useMaxY) My = options.maxY;

      // Set the bounds in the graph options object
      options.hasBounds = true;
      options.minX = mx;
      options.maxX = Mx;
      options.minY = my;
      options.maxY = My;
    }
    else {
      mx = options.minX;
      Mx = options.maxX;
      my = options.minY;
      My = options.maxY;
    }
    if (options.useBcgd) background = options.bcgd;

    // Set widths
    float wx = Mx - mx;
    float wy = My - my;

    // Set the background
    if (options.paintBackground) coverPalette(background);
    options.paintBackground = false;

    auto colorF = [=] (float, float, bool&) -> RGBApixel {
      return lineColor;
    };
      
    // Draw graph
    pair<float, float> last = first;
    auto p = data.begin();
    float xf0, xf1, yf0, yf1;
    xf0 = (p->first  - mx)/wx;
    yf0 = (p->second - my)/wy;
    ++p; // Increment p
    // Go through all data
    for (; p!=data.end(); ++p) {
      float xf1 = (p->first - mx)/wx;
      float yf1 = (p->second - my)/wy;
      drawLineByFactors(xf0, yf0, xf1, yf1, colorF, 0);
      xf0 = xf1;
      yf0 = yf1;
    }
  }

  void Palette::pixelsToImage(BMP *photo, int mx, int Mx, int my, int My) {
    if (photo==nullptr) return;
    // Set photo's pixels
    for (int y=my; y<My; ++y)
      for (int x=mx; x<Mx; ++x) {
        // Get the pixel data 
        PixL &p = pixelData->at(x, y);
        if (p.empty()) photo->SetPixel(x-mx, y-my, RGB_Black);
        else {
          // Combine pixel data
          switch (combinationRule) {
            default:
            case 0: {
              photo->SetPixel(x-mx, y-my, p.first());
              break;
            }
            case 1: {
              int r(0), g(0), b(0);
              for (auto c : p.data) {
                if (c.second.Red>r)   r = c.second.Red;
                if (c.second.Green>g) g = c.second.Green;
                if (c.second.Blue>b)  b = c.second.Blue;
              }
              photo->SetPixel(x-mx, y-my, RGBApixel(r, g, b));
              break;
            }
            case 2: {
              photo->SetPixel(x-mx, y-my, RGBApixel(p.size()*80, 155, p.size()*100));
              break;
            }
          }
        }
      }
  }

}