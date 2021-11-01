"use strict";

const Filters = {};

////////////////////////////////////////////////////////////////////////////////
// General utility functions
////////////////////////////////////////////////////////////////////////////////

// Hardcoded Pi value
// const pi = 3.14159265359;
const pi = Math.PI;

// Constrain val to the range [min, max]
function clamp(val, min, max) {
    /* Shorthand for:
    * if (val < min) {
    *   return min;
    * } else if (val > max) {
    *   return max;
    * } else {
    *   return val;
    * }
    */
    return val < min ? min : val > max ? max : val;
}

// Extract vertex coordinates from a URL string
function stringToCoords(vertsString) {
    const centers = [];
    const coordStrings = vertsString.split("x");
    for (let i = 0; i < coordStrings.length; i++) {
        const coords = coordStrings[i].split("y");
        const x = parseInt(coords[0]);
        const y = parseInt(coords[1]);
        if (!isNaN(x) && !isNaN(y)) {
            centers.push({ x: x, y: y });
        }
    }

    return centers;
}

// Blend scalar start with scalar end. Note that for image blending,
// end would be the upper layer, and start would be the background
function blend(start, end, alpha) {
    return start * (1 - alpha) + end * alpha;
}

// ----------- STUDENT CODE BEGIN ------------
// ----------- Our reference solution uses 72 lines of code.
// ----------- STUDENT CODE END ------------

////////////////////////////////////////////////////////////////////////////////
// Filters
////////////////////////////////////////////////////////////////////////////////

// You've already implemented this in A0! Feel free to copy your code into here
Filters.fillFilter = function(image, color) {
    image.fill(color);

    return image;
};

// You've already implemented this in A0! Feel free to copy your code into here
Filters.brushFilter = function(image, radius, color, vertsString) {
    // centers is an array of (x, y) coordinates that each defines a circle center
    const centers = stringToCoords(vertsString);

    // draw a filled circle centered at every location in centers[].
    // radius and color are specified in function arguments.
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 10 lines of code.
    for (var i = 0; i < centers.length; i++) {
        for (var x = centers[i].x - radius; x < centers[i].x + radius; x++) {
            for (var y = centers[i].y - radius; y < centers[i].y + radius; y++) {
                var in_circle = (y - centers[i].y) ** 2 + (x - centers[i].x) ** 2;
                if (in_circle <= radius ** 2) {
                    image.setPixel(x, y, color);
                }
            }
        }
    }
    // ----------- STUDENT CODE END ------------

    return image;
};

// You've already implemented this in A0! Feel free to copy your code into here
Filters.softBrushFilter = function(image, radius, color, alpha_at_center, vertsString) {
    // centers is an array of (x, y) coordinates that each defines a circle center
    const centers = stringToCoords(vertsString);

    // draw a filled circle with opacity equals to alpha_at_center at the center of each circle
    // the opacity decreases linearly along the radius and becomes zero at the edge of the circle
    // radius and color are specified in function arguments.
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 20 lines of code.
    var original_color;
    var factor;
    for (var i = 0; i < centers.length; i++) {
        for (var x = centers[i].x - radius; x < centers[i].x + radius; x++) {
            for (var y = centers[i].y - radius; y < centers[i].y + radius; y++) {
                var in_circle = (y - centers[i].y) ** 2 + (x - centers[i].x) ** 2;
                if (in_circle <= radius ** 2) {
                    factor = (radius - in_circle ** (1 / 2)) / radius
                    original_color = image.getPixel(x, y);
                    original_color.data[0] = (1 - alpha_at_center * factor) * original_color.data[0] + alpha_at_center * color.data[0] * factor;
                    original_color.data[1] = (1 - alpha_at_center * factor) * original_color.data[1] + alpha_at_center * color.data[1] * factor;
                    original_color.data[2] = (1 - alpha_at_center * factor) * original_color.data[2] + alpha_at_center * color.data[2] * factor;
                    image.setPixel(x, y, original_color);
                }
            }
        }
    }
    // ----------- STUDENT CODE END ------------

    return image;
};

// Ratio is a value in the domain [-1, 1]. When ratio is < 0, linearly blend the image
// with black. When ratio is > 0, linearly blend the image with white. At the extremes
// of -1 and 1, the image should be completely black and completely white, respectively.
Filters.brightnessFilter = function(image, ratio) {
    let alpha, dirLuminance;
    if (ratio < 0.0) {
        alpha = 1 + ratio;
        dirLuminance = 0; // blend with black
    } else {
        alpha = 1 - ratio;
        dirLuminance = 1; // blend with white
    }

    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);

            pixel.data[0] = clamp(alpha * pixel.data[0] + (1 - alpha) * dirLuminance, 0, 1);
            pixel.data[1] = clamp(alpha * pixel.data[1] + (1 - alpha) * dirLuminance, 0, 1);
            pixel.data[2] = clamp(alpha * pixel.data[2] + (1 - alpha) * dirLuminance, 0, 1);

            image.setPixel(x, y, pixel);
        }
    }

    return image;
};

// Reference at this:
//      https://en.wikipedia.org/wiki/Image_editing#Contrast_change_and_brightening
// value = (value - 0.5) * (tan ((contrast + 1) * PI/4) ) + 0.5;
// Note that ratio is in the domain [-1, 1]
Filters.contrastFilter = function(image, ratio) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 14 lines of code.
    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);

            pixel.data[0] = (pixel.data[0] - 0.5) * (Math.tan((ratio + 1) * Math.PI/4)) + 0.5;
            pixel.data[1] = (pixel.data[1] - 0.5) * (Math.tan((ratio + 1) * Math.PI/4)) + 0.5;
            pixel.data[2] = (pixel.data[2] - 0.5) * (Math.tan((ratio + 1) * Math.PI/4)) + 0.5;

            image.setPixel(x, y, pixel);
        }
    }
    // ----------- STUDENT CODE END ------------
    return image;
};

// Note that the argument here is log(gamma)
Filters.gammaFilter = function(image, logOfGamma) {
    const gamma = Math.exp(logOfGamma);
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 9 lines of code.
    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);

            pixel.data[0] = pixel.data[0] ** gamma;
            pixel.data[1] = pixel.data[1] ** gamma;
            pixel.data[2] = pixel.data[2] ** gamma;

            image.setPixel(x, y, pixel);
        }
    }
    // ----------- STUDENT CODE END ------------
    //Gui.alertOnce ('gammaFilter is not implemented yet');
    return image;
};

/*
* The image should be perfectly clear up to innerRadius, perfectly dark
* (black) at outerRadius and beyond, and smoothly increase darkness in the
* circular ring in between. Both are specified as multiples of half the length
* of the image diagonal (so 1.0 is the distance from the image center to the
* corner).
*
* Note that the vignette should still form a perfect circle!
*/
Filters.vignetteFilter = function(image, innerR, outerR) {
    // Let's ensure that innerR is at least 0.1 smaller than outerR
    innerR = clamp(innerR, 0, outerR - 0.1);
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 17 lines of code.
    var halfDiag = ((image.width**2 + image.height**2)**.5)*.5;
    var xCenter = image.width/2;
    var yCenter = image.height/2;
    var distance;
    var angle;

    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);
            distance = (((x-xCenter)**2 + (y-yCenter)**2)**.5)/halfDiag; // distance from center normalized by half diagonal

            // if in inner radius, no changes
            if(distance < innerR){
                continue
            }

            // if between inner and outer radius, apply the natural vignetting filter
            else if(distance > innerR && distance < outerR) {
                let n = 1-(distance-innerR)/(outerR - innerR);
                pixel.data[0] = pixel.data[0]*n;
                pixel.data[1] = pixel.data[1]*n;
                pixel.data[2] = pixel.data[2]*n;
            }

            // if outside outer radius, set pixel value to 0
            else{
                pixel.data[0] = 0;
                pixel.data[1] = 0;
                pixel.data[2] = 0;
            }

            image.setPixel(x, y, pixel);
        }
    }
    // ----------- STUDENT CODE END ------------
    //Gui.alertOnce ('vignetteFilter is not implemented yet');
    return image;
};

/*
* You will want to build a normalized CDF of the L channel in the image.
*/
Filters.histogramEqualizationFilter = function(image) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 33 lines of code.



    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('histogramEqualizationFilter is not implemented yet');
    return image;
};

// Set each pixel in the image to its luminance
Filters.grayscaleFilter = function(image) {
    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);
            const luminance = 0.2126 * pixel.data[0] + 0.7152 * pixel.data[1] + 0.0722 * pixel.data[2];
            pixel.data[0] = luminance;
            pixel.data[1] = luminance;
            pixel.data[2] = luminance;

            image.setPixel(x, y, pixel);
        }
    }

    return image;
};

// Adjust each channel in each pixel by a fraction of its distance from the average
// value of the pixel (luminance).
// See: http://www.graficaobscura.com/interp/index.html
Filters.saturationFilter = function(image, ratio) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 13 lines of code.
    const greyImage = this.grayscaleFilter(image.copyImg());
    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);
            const greyPixel = greyImage.getPixel(x,y);

            pixel.data[0] = pixel.data[0]  + (pixel.data[0] - greyPixel.data[0]) * ratio;
            pixel.data[1] = pixel.data[1]  + (pixel.data[1] - greyPixel.data[1]) * ratio;
            pixel.data[2] = pixel.data[2]  + (pixel.data[2] - greyPixel.data[2]) * ratio;
            pixel.clamp()

            image.setPixel(x, y, pixel);
        }
    }
    // ----------- STUDENT CODE END ------------
    //Gui.alertOnce ('saturationFilter is not implemented yet');
    return image;
};

// Apply the Von Kries method: convert the image from RGB to LMS, divide by
// the LMS coordinates of the white point color, and convert back to RGB.
Filters.whiteBalanceFilter = function(image, white) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 23 lines of code.

    // convert reference white to lms colorspace
    white = white.rgbToXyz();
    white = white.xyzToLms();

    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            // convert pixel to lms colorspace
            var pixel = image.getPixel(x, y);
            pixel = pixel.rgbToXyz();
            pixel = pixel.xyzToLms();
            
            // apply white balance
            pixel.data[0] = pixel.data[0]/white.data[0];
            pixel.data[1] = pixel.data[1]/white.data[1];
            pixel.data[2] = pixel.data[2]/white.data[2];

            // convert back to rgb colorspace
            pixel = pixel.lmsToXyz();
            pixel = pixel.xyzToRgb();
            
            // clamp pixel values
            pixel.clamp();
            
            image.setPixel(x, y, pixel);
        }
    }
    // ----------- STUDENT CODE END ------------
    return image;
};

// This is similar to the histogram filter, except here you should take the
// the CDF of the L channel in one image and
// map it to another
//
Filters.histogramMatchFilter = function(image, refImg) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 58 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('histogramMatchFilter is not implemented yet');
    return image;
};

// Convolve the image with a gaussian filter.
// NB: Implement this as a seperable gaussian filter
Filters.gaussianFilter = function(image, sigma) {
    // note: this function needs to work in a new copy of the image
    //       to avoid overwriting original pixels values needed later
    // create a new image with the same size as the input image
    let newImg = image.createImg(image.width, image.height);
    // the filter window will be [-winR, winR] for a total diameter of roughly Math.round(3*sigma)*2+1;
    const winR = Math.round(sigma * 3);
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 58 lines of code.
    
    // initialize 
    const kernel = [];
    for (let i = -1 * winR; i <= winR; i++){
        kernel.push(Math.exp(-1 * (i ** 2) / (2*sigma)) / ((2*Math.PI * sigma**2)**.5));
    }

    let normal = kernel.reduce((a,b) => {return a + b});

    for (let i = 0; i < kernel.length; i++){
        kernel[i] = kernel[i]/ normal;
    }
    // Convolve horizontally
    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            var rSum, gSum, bSum;
            rSum = gSum = bSum = 0;
            for (let i = 0; i < kernel.length; i++){
                let pixel = image.getPixel(clamp(x - winR + i, 0, image.width), y);
                rSum += pixel.data[0] * kernel[i];
                gSum += pixel.data[1] * kernel[i];
                bSum += pixel.data[2] * kernel[i];
            }
            let new_pixel = new Pixel(rSum, gSum, bSum);
            newImg.setPixel(x, y, new_pixel);
        }
    }

    // Convolve vertically
    for (let x = 0; x < newImg.width; x++) {
        for (let y = 0; y < newImg.height; y++) {
            let rSum, gSum, bSum;
            rSum = gSum = bSum = 0;
            for (let i = 0; i < kernel.length; i++){
                let pixel = newImg.getPixel(x, clamp(y - winR + i, 0, newImg.height));
                rSum += pixel.data[0] * kernel[i];
                gSum += pixel.data[1] * kernel[i];
                bSum += pixel.data[2] * kernel[i];
            }
            let new_pixel = new Pixel(rSum, gSum, bSum);
            newImg.setPixel(x, y, new_pixel);
        }
    }

    // ----------- STUDENT CODE END ------------
    //Gui.alertOnce ('gaussianFilter is not implemented yet');
    return newImg;
};

/*
* First the image with the edge kernel and then add the result back onto the
* original image.
*/
Filters.sharpenFilter = function(image) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 33 lines of code.
    const kernel = [];
    const winSize = 1;
    let newImg = image.copyImg();

    // sharpening kernel intialization
    for(let i = 0; i < 3; i++){
        kernel.push([])
        for(let j = 0; j < 3; j++){
            kernel[i].push(-1);
        }
    }
    kernel[1][1] = 9;

    // convolve sharpen filter
    for (let x = 0; x < newImg.width; x++) {
        for (let y = 0; y < newImg.height; y++) {
            let rSum, gSum, bSum;
            rSum = gSum = bSum = 0;
            for (let i = 0; i < kernel.length; i++){
                for (let j = 0; j < kernel.length; j++){
                    let pixel = newImg.getPixel(clamp(x - winSize + i, 0, newImg.width), clamp(y - winSize + j, 0, newImg.height));
                    rSum += pixel.data[0] * kernel[i][j];
                    gSum += pixel.data[1] * kernel[i][j];
                    bSum += pixel.data[2] * kernel[i][j];
                }
            }
            let new_pixel = new Pixel(rSum, gSum, bSum);
            image.setPixel(x, y, new_pixel);
        }
    }

    // ----------- STUDENT CODE END ------------
    return image;
};

/*
* Convolve the image with the edge kernel from class. You might want to define
* a convolution utility that convolves an image with some arbitrary input kernel
*
* For this filter, we recommend inverting pixel values to enhance edge visualization
*/
Filters.edgeFilter = function(image) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 57 lines of code.
    const kernel = [];
    const winSize = 1;
    let newImg = image.copyImg();

    // sharpening kernel intialization
    for(let i = 0; i < 3; i++){
        kernel.push([])
        for(let j = 0; j < 3; j++){
            kernel[i].push(-1);
        }
    }
    kernel[1][1] = 8;

    // convolve sharpen filter
    for (let x = 0; x < newImg.width; x++) {
        for (let y = 0; y < newImg.height; y++) {
            let rSum, gSum, bSum;
            rSum = gSum = bSum = 0;
            for (let i = 0; i < kernel.length; i++){
                for (let j = 0; j < kernel.length; j++){
                    let pixel = newImg.getPixel(clamp(x - winSize + i, 0, newImg.width), clamp(y - winSize + j, 0, newImg.height));
                    rSum += pixel.data[0] * kernel[i][j];
                    gSum += pixel.data[1] * kernel[i][j];
                    bSum += pixel.data[2] * kernel[i][j];
                }
            }
            // implementing 1 - pixel for visualization
            let new_pixel = new Pixel(1 - rSum, 1 - gSum, 1 - bSum);
            image.setPixel(x, y, new_pixel);
        }
    }
    // ----------- STUDENT CODE END ------------
    return image;
};

// Set a pixel to the median value in its local neighbor hood. You might want to
// apply this seperately to each channel.
Filters.medianFilter = function(image, winR) {
    // winR: the window will be  [-winR, winR];
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 36 lines of code.
    var localR, localG, localB;
    localR = localB = localG = [];
    let newImg = image.copyImg();

    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            localR = [];
            localB = [];
            localG = [];

            for (let i = 0; i < 2 * winR + 1; i++){
                for (let j = 0; j < 2 * winR + 1; j++){
                    let pixel = newImg.getPixel(clamp(x - winR + i, 0, newImg.width), clamp(y - winR + j, 0, newImg.height));
                    localR.push(pixel.data[0]);
                    localG.push(pixel.data[1]);
                    localB.push(pixel.data[2]);
                }
            }
            // sort and choose median pixel
            let new_pixel = new Pixel(localR.sort(function(a, b){return b - a})[winR], localG.sort(function(a, b){return b - a})[winR], localB.sort(function(a, b){return b - a})[winR]);
            image.setPixel(x, y, new_pixel);
        }
    }

    // ----------- STUDENT CODE END ------------
    return image;
};

// Apply a bilateral filter to the image. You will likely want to reference
// precept slides, lecture slides, and the assignments/examples page for help.
Filters.bilateralFilter = function(image, sigmaR, sigmaS) {
    // reference: https://en.wikipedia.org/wiki/Bilateral_filter
    // we first compute window size and preprocess sigmaR
    const winR = Math.round((sigmaR + sigmaS) * 1.5);
    sigmaR = sigmaR * (Math.sqrt(2) * winR);


    let newImg = image.copyImg();
    const window = Math.round(2*Math.max(sigmaR,sigmaS)) *2 +1;

    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 53 lines of code.

    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            let localR = 0;
            let localB = 0;
            let localG = 0;
            let rNormal = 0;
            let gNormal = 0;
            let bNormal = 0;
            let pixel = newImg.getPixel(x,y);
            for (let i = 0; i < 2*winR+1; i++){
                for (let j = 0; j < 2*winR+1; j++){
                    let xLocal = clamp(x - winR + i, 0, newImg.width)
                    let yLocal = clamp(y - winR + j, 0, newImg.height);
                    let localPixel = newImg.getPixel(xLocal, yLocal);
                    let spacialFactor = ((x-xLocal)**2 + (y-yLocal)**2)/(2*(sigmaS**2));
                    
                    // red color
                    let colorFactor = ((pixel.data[0] - localPixel.data[0])**2)/(2*(sigmaR**2));
                    localR += (localPixel.data[0]*Math.exp(-1*(spacialFactor + colorFactor)));
                    rNormal += Math.exp(-1*(spacialFactor + colorFactor))

                    // green color
                    colorFactor = ((pixel.data[1] - localPixel.data[1])**2)/(2*(sigmaR**2));
                    localG += (localPixel.data[1]*Math.exp(-1*(spacialFactor + colorFactor)));
                    gNormal += Math.exp(-1*(spacialFactor + colorFactor))

                    //blue color
                    colorFactor = ((pixel.data[2] - localPixel.data[2])**2)/(2*(sigmaR**2));
                    localB += (localPixel.data[2]*Math.exp(-1*(spacialFactor + colorFactor)));
                    bNormal += Math.exp(-1*(spacialFactor + colorFactor))
                }
            }
            //let sum = localR.reduce((a,b) => {return a + b});
            let r = localR/rNormal;
            let g = localG/gNormal;
            let b = localB/bNormal;

            // sort and choose median pixel
            let new_pixel = new Pixel(r, g, b);
            image.setPixel(x, y, new_pixel);
        }
    }

    // ----------- STUDENT CODE END ------------
    //Gui.alertOnce ('bilateralFilter is not implemented yet');
    return image;
};

// Conver the image to binary
Filters.quantizeFilter = function(image) {
    // convert to grayscale
    image = Filters.grayscaleFilter(image);

    // use center color
    for (let i = 0; i < image.height; i++) {
        for (let j = 0; j < image.width; j++) {
            const pixel = image.getPixel(j, i);
            for (let c = 0; c < 3; c++) {
                pixel.data[c] = Math.round(pixel.data[c]);
            }
            pixel.clamp();
            image.setPixel(j, i, pixel);
        }
    }
    return image;
};

// To apply random dithering, first convert the image to grayscale, then apply
// random noise, and finally quantize
Filters.randomFilter = function(image) {
    // convert to grayscale
    image = Filters.grayscaleFilter(image);

    for (let i = 0; i < image.height; i++) {
        for (let j = 0; j < image.width; j++) {
            const pixel = image.getPixel(j, i);
            let n = Math.random() - .5;
            for (let c = 0; c < 3; c++) {
                pixel.data[c] = Math.round(pixel.data[c] - n);
            }
            pixel.clamp();
            image.setPixel(j, i, pixel);
        }
    } 

    return image;
};

// Apply the Floyd-Steinberg dither with error diffusion
Filters.floydFilter = function(image) {
    // convert to grayscale
    image = Filters.grayscaleFilter(image);

    for (let i = 0; i < image.height; i++) {
        for (let j = 0; j < image.width; j++) {
            const pixel = image.getPixel(j, i);
            
            // Retrieve neighboring pixels
            const pixel_r = image.getPixel(Math.min(j + 1, image.width - 1), i);
            const pixel_br = image.getPixel(Math.min(j + 1, image.width - 1), Math.min(i + 1, image.height -1));
            const pixel_b = image.getPixel(j, Math.min(i + 1, image.height - 1));
            const pixel_bl = image.getPixel(Math.max(j - 1, 0), Math.min(i + 1, image.height -1));
            
            for (let c = 0; c < 3; c++) {
                var p = Math.round(pixel.data[c]);
                var q = pixel.data[c] - p;
                pixel.data[c] = p;
                
                // Diffuse error
                pixel_r.data[c] += 7*q/16;
                pixel_br.data[c] += q/16;
                pixel_b.data[c] += 5*q/16;
                pixel_bl.data[c] += 3*q/16;
            }
            pixel.clamp();
            image.setPixel(j, i, pixel);

            // Update neighboring pixels
            image.setPixel(Math.min(j + 1, image.width - 1), i, pixel_r);
            image.setPixel(Math.min(j + 1, image.width - 1), Math.min(i + 1, image.height - 1), pixel_br);
            image.setPixel(j, Math.min(i + 1, image.height - 1), pixel_b);
            image.setPixel(Math.max(j - 1, 0), Math.min(i + 1, image.height - 1), pixel_bl);
        }
    }

    return image;
};

// Apply ordered dithering to the image. We recommend using the pattern from the
// examples page and precept slides.
Filters.orderedFilter = function(image) {
    // convert to gray scale
    image = Filters.grayscaleFilter(image);
    
    // Bayer 4 pattern 
    const m = [
        [15, 7, 13, 5],
        [3, 11, 1, 9],
        [12, 4, 14, 6], 
        [0, 8, 2, 10]
    ];

    for (let i = 0; i < image.height; i++) {
        for (let j = 0; j < image.width; j++) {
            const pixel = image.getPixel(j, i);
            var x = j % 4;
            var y = i % 4;
            var threshold = (m[x][y] + 1) / (m.length**2 + 1);
            for (let c = 0; c < 3; c++) {
                var err = pixel.data[c];
                if (err < threshold) pixel.data[c] = 0;
                else pixel.data[c] = 1;
            }
            pixel.clamp();
            image.setPixel(j, i, pixel);
        }
    } 
    
    return image;
};

// Implement bilinear and Gaussian sampling (in addition to the basic point sampling).
// This operation doesn't appear on GUI and should be used as a utility function.
// Call this function from filters that require sampling (e.g. scale, rotate)
Filters.samplePixel = function(image, x, y, mode) {

    if (mode === "bilinear") {
        if (x >= image.width || y >= image.width || x < 0 || y < 0){
            y = Math.max(0, Math.min(Math.round(y), image.height - 1));
            x = Math.max(0, Math.min(Math.round(x), image.width - 1));
            return image.getPixel(x, y);
        }
        const pixel = new Pixel(0, 0, 0);

        // Find nearest integer coordinates and retrieve pixels
        var x1 = Math.floor(x);
        var x2 = Math.ceil(x);
        var y1 = Math.floor(y); 
        var y2 = Math.ceil(y);

        if (x1 == x && y1 == y) return image.getPixel(x, y);   
        else if (x1 == x) {
            const p1 = image.getPixel(x, y1);
            const p2 = image.getPixel(x, y2);
            for (let c = 0; c < 3; c++) {
                pixel.data[c] = (p1.data[c]*(y2-y)) + p2.data[c]*(y-y1);
            }
            return pixel;
        } 
        else if (y1 == y) {
            const p1 = image.getPixel(x1, y);
            const p2 = image.getPixel(x2, y);
            for (let c = 0; c < 3; c++) { 
                pixel.data[c] = (p1.data[c]*(x2-x)) + p2.data[c]*(x-x1);
            }
            return pixel;
        }                 

        const p11 = image.getPixel(x1, y1);
        const p12 = image.getPixel(x1, y2);
        const p21 = image.getPixel(x2, y1);
        const p22 = image.getPixel(x2, y2);

        for (let c = 0; c < 3; c++) {
            pixel.data[c] = (p11.data[c]*(x2-x)*(y2-y) + p21.data[c]*(x-x1)*(y2-y) + p12.data[c]*(x2-x)*(y-y1) + p22.data[c]*(x-x1)*(y-y1));
        }

        return pixel;

    } else if (mode === "gaussian") {
        // ----------- Our reference solution uses 38 lines of code.
        if (x >= image.width || y >= image.width || x < 0 || y < 0){
            y = Math.max(0, Math.min(Math.round(y), image.height - 1));
            x = Math.max(0, Math.min(Math.round(x), image.width - 1));
            return image.getPixel(x, y);
        }
        const pixel = new Pixel(0, 0, 0);
        const x1 = Math.round(x);
        const y1 = Math.round(y);
        const sigma = 1;
        const normalization = [0,0,0]

        for (let x2 = -3*sigma; x2 < 3*sigma + 1; x2++) {
            for (let y2 = -3*sigma; y2 < 3*sigma + 1; y2++) {

                var cur_x = clamp(x2 + x1, 0, image.width+1);
                var cur_y = clamp(y2 + y1, 0, image.height+1);

                const curPixel = image.getPixel(cur_x, cur_y);
                var d = (x-cur_x)**2 + (y-cur_y)**2;

                for (let c = 0; c < 3; c++) {
                    pixel.data[c] += curPixel.data[c]*Math.exp((-(d)/(2*(sigma**2))));
                    normalization[c] += Math.exp((-(d)/(2*(sigma**2))));
                }

                pixel.data[0] /= normalization[0];
                pixel.data[1] /= normalization[1];
                pixel.data[2] /= normalization[2];
                return pixel;
            }
        }
        

    } else {
        // point sampling
        y = Math.max(0, Math.min(Math.round(y), image.height - 1));
        x = Math.max(0, Math.min(Math.round(x), image.width - 1));
        return image.getPixel(x, y);
    }
};

// Translate the image by some x, y and using a requested method of sampling/resampling
Filters.translateFilter = function(image, x, y, sampleMode) {
    // Note: set pixels outside the image to RGBA(0,0,0,0)
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 21 lines of code.
    const newImg = image.copyImg();
 
    for (let x1 = 0; x1 < image.width; x1++) {
        for (let y1 = 0; y1 < image.height; y1++) {
            let dx =  x1 - x;
            let dy = y1 - y;
            if (dx < 0 || dy < 0 || dx > image.width || dy > image.height){
                let pixel = new Pixel(0,0,0);
                newImg.setPixel(x1,y1,pixel);
            }
            else {
                const pixel = this.samplePixel(image, x1 - x, y1 - y, sampleMode); 
                newImg.setPixel(x1, y1, pixel);
            }
        }
    } 

    // ----------- STUDENT CODE END ------------
    return newImg;
};

// Scale the image by some ratio and using a requested method of sampling/resampling
Filters.scaleFilter = function(image, ratio, sampleMode) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 19 lines of code.
    const newImg = image.createImg(Math.round(image.width*ratio), Math.round(image.height*ratio));
    
    for (let x1 = 0; x1 < newImg.width; x1++) {
        for (let y1 = 0; y1 < newImg.height; y1++) {
            const pixel = this.samplePixel(image, x1/ratio, y1/ratio, sampleMode); 
            newImg.setPixel(x1, y1, pixel);
        }
    } 


    // ----------- STUDENT CODE END ------------
   // Gui.alertOnce ('scaleFilter is not implemented yet');
    return newImg;
};

// Rotate the image by some angle and using a requested method of sampling/resampling
Filters.rotateFilter = function(image, radians, sampleMode) {
    // Note: set pixels outside the image to RGBA(0,0,0,0)
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 29 lines of code.
    // ----------- STUDENT CODE END ------------

    let new_width = 0;
    let new_height = 0;
    if ((0 <= radians && radians < Math.PI/2) || (Math.PI < radians && radians < 3*Math.PI/2)) {
        new_width = Math.round(image.width*Math.cos(radians) + image.height*Math.sin(radians));
        new_height = Math.round(image.width*Math.sin(radians) + image.height*Math.cos(radians));
    }
    else {
        let w = image.height;
        let h = image.width; 
        let rad = radians - (Math.PI/2);

        new_width = Math.round(w*Math.cos(rad) + h*Math.sin(rad));
        new_height = Math.round(w*Math.sin(rad) + h*Math.cos(rad));
    }
    
    const newImg = image.createImg(new_width, new_height);

    for (let x = 0; x < new_width; x++) {
        for (let y = 0; y < new_height; y++) {
            let x_difference = (x + (new_width - image.width)) * Math.cos(radians) - (y)*Math.sin(radians); // +73
            let y_difference = (x) * Math.sin(radians) + (y-(new_height - image.height))*Math.cos(radians) ;
            if (x_difference > image.width || y_difference > image.height || x_difference < 0|| y_difference < 0){
                let pixel = new Pixel(0,0,0);
                newImg.setPixel(x, y, pixel);
            }
            else{
                let pixel = this.samplePixel(image,  x_difference, y_difference, sampleMode);
                newImg.setPixel(x, y, pixel); 
            }
            //newImg.setPixel(x, y, pixel);
        }
    } 

    return newImg;
};

// Swirl the filter about its center. The rotation of the swirl should be in linear increase
// along the radial axis up to radians
Filters.swirlFilter = function(image, radians, sampleMode) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 26 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('swirlFilter is not implemented yet');
    return image;
};

// Set alpha from luminance
Filters.getAlphaFilter = function(backgroundImg, foregroundImg) {
    for (let i = 0; i < backgroundImg.height; i++) {
        for (let j = 0; j < backgroundImg.width; j++) {
            const pixelBg = backgroundImg.getPixel(j, i);
            const pixelFg = foregroundImg.getPixel(j, i);
            const luminance =
            0.2126 * pixelFg.data[0] + 0.7152 * pixelFg.data[1] + 0.0722 * pixelFg.data[2];
            pixelBg.a = luminance;
            backgroundImg.setPixel(j, i, pixelBg);
        }
    }

    return backgroundImg;
};

// Composites the foreground image over the background image, using the alpha
// channel of the foreground image to blend two images.
Filters.compositeFilter = function(backgroundImg, foregroundImg) {
    // Assume the input images are of the same sizes.
    for (let i = 0; i < backgroundImg.height; i++) {
        for (let j = 0; j < backgroundImg.width; j++) {
            const pixelBg = backgroundImg.getPixel(j, i);
            const pixelFg = foregroundImg.getPixel(j, i);
            const alpha = pixelFg.a;
            for (let c = 0; c < 3; c++) {
                pixelBg.data[c] = alpha*pixelFg.data[c] + (1-alpha)*pixelBg.data[c];
            }
            backgroundImg.setPixel(j, i, pixelBg);
        }
    }

    return backgroundImg;
};

Filters.warpPixel = function(lines1, lines2, x, y) {
    
    const p = .5;
    const a = .01;
    const b = 2;

    let dsum = [0,0];
    let weightsum = 0;
    
    // Iterate through all lines
    for (let i = 0; i < lines1.length; i++){
        let dx = lines1[i].x1 - lines1[i].x0;
        let dy = lines1[i].y1 - lines1[i].y0;
        let px = x - lines1[i].x0;
        let py = y - lines1[i].y0;
        let perpy = -(lines1[i].x1 - lines1[i].x0);
        let perpx = lines1[i].y1 - lines1[i].y0;

        let u = (px*dx + py*dy) / (dx**2 + dy**2);
        let v = (px*perpx + py*perpy) / Math.sqrt(dx**2 + dy**2);
        
        let perpx_prime = lines2[i].y1 - lines2[i].y0;
        let perpy_prime = -(lines2[i].x1 - lines2[i].x0);
        let dx_prime = lines2[i].x1 - lines2[i].x0;
        let dy_prime = lines2[i].y1 - lines2[i].y0;

        let x_prime = lines2[i].x0 + u * (dx_prime) + v * perpx_prime/Math.sqrt(dx_prime**2 + dy_prime**2);
        let y_prime = lines2[i].y0 + u * (dy_prime) + v * perpy_prime/Math.sqrt(dy_prime**2 + dx_prime**2);

        let displacement = [(x_prime - x), (y_prime - y)];
        let dist;

        if(u >= 0 && u <= 1){
            dist = Math.abs(v);
        }
        else if(u < 0){
            dist = Math.sqrt(px**2 + py**2);
        }
        else if(u > 1){
            let qx = x - lines1[i].x1;
            let qy = y - lines1[i].y1;
            dist = Math.sqrt(qx**2 + qy**2);
        }

        let length = Math.sqrt(dx_prime**2 + dy_prime**2);
        let weight = ((length ** p)/ (a+dist))**b;
        
        dsum[0] += displacement[0]*weight;
        dsum[1] += displacement[1]*weight;
        weightsum += weight;
    }
    let X = x + dsum[0]/weightsum;
    let Y = y + dsum[1]/weightsum;

    return [X, Y];
}

// Morph two images according to a set of correspondance lines
Filters.morphFilter = function(initialImg, finalImg, alpha, sampleMode, linesFile) {
    const lines = Parser.parseJson("images/" + linesFile);

    // The provided linesFile represents lines in a flipped x, y coordinate system
    //  (i.e. x for vertical direction, y for horizontal direction).
    // Therefore we first fix the flipped x, y coordinates here.
    for (let i = 0; i < lines.initial.length; i++) {
        [lines.initial[i].x0, lines.initial[i].y0] = [lines.initial[i].y0, lines.initial[i].x0];
        [lines.initial[i].x1, lines.initial[i].y1] = [lines.initial[i].y1, lines.initial[i].x1];
        [lines.final[i].x0, lines.final[i].y0] = [lines.final[i].y0, lines.final[i].x0];
        [lines.final[i].x1, lines.final[i].y1] = [lines.final[i].y1, lines.final[i].x1];
    }
    
    // ----------- Our reference solution uses 114 lines of code.
    const image = finalImg.copyImg();
    var current_lines = [];
    
    for(let i = 0; i < lines.initial.length; i++){
        const x0 = alpha*lines.initial[i].x0 + (1-alpha)*lines.final[i].x0
        const y0 = alpha*lines.initial[i].y0 + (1-alpha)*lines.final[i].y0;
        const x1 = alpha*lines.initial[i].x1 + (1-alpha)*lines.final[i].x1;
        const y1 = alpha*lines.initial[i].y1 + (1-alpha)*lines.final[i].y1;
        current_lines.push({x0: x0, y0: y0, x1: x1, y1: y1});
    }

    // Iterate through all pixels
    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++){
            let warp1 = this.warpPixel(lines.initial, current_lines, x, y);
            let warp2 = this.warpPixel(lines.final, current_lines, x, y);

            let initial_pixel = this.samplePixel(initialImg, warp2[0], warp2[1], sampleMode); 
            let final_pixel = this.samplePixel(finalImg, warp1[0], warp1[1], sampleMode); 
            
            let r = (1-alpha) * initial_pixel.data[0] + alpha * final_pixel.data[0];
            let g = (1-alpha) * initial_pixel.data[1] + alpha * final_pixel.data[1];
            let b = (1-alpha) * initial_pixel.data[2] + alpha * final_pixel.data[2];

            let new_pixel = new Pixel(r,g,b);
            image.setPixel(x, y, new_pixel)
        }
    } 

    return image;
};

// Use k-means to extract a pallete from an image
Filters.paletteFilter = function(image, colorNum) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 89 lines of code.
    // initialize centroids for k-means clustering
    const centroids = [];
    let centroid_elements = [];
    let basecolors = [new Pixel(1,1,1), new Pixel(0,1,0), new Pixel(0, 0, 1), new Pixel(1,0,0), new Pixel(0,1,1), new Pixel(1,0,1)];
    for (let i = 0; i < colorNum; i++){
        centroids.push(basecolors[i]);
        centroid_elements.push([]);
    }

    // iterate through clustering algorithm
    for(let iterations = 0; iterations < 10; iterations++){
        centroid_elements = [];
        for (let i = 0; i < colorNum; i++){
            centroid_elements.push([]);
        }
        // Look at every pixel in image
        for (let x = 0; x < image.width; x++) {
            for (let y = 0; y < image.height; y++){
                let p = image.getPixel(x,y);
                let best_dist = (centroids[0].data[0] - p.data[0])**2 +(centroids[0].data[1] - p.data[1])**2 +(centroids[0].data[2] - p.data[2])**2;
                let best_point = 0

                // find which centroid the pixel is most closely associated with
                for(let k = 0; k < centroids.length; k++){
                    let dist = (centroids[k].data[0] - p.data[0])**2 +(centroids[k].data[1] - p.data[1])**2 +(centroids[k].data[2] - p.data[2])**2;
                    if(dist < best_dist){
                        best_dist = dist;
                        best_point = k;
                    }
                }

                centroid_elements[best_point].push(p);
            }
        } 

        // Recalulate centers
        for (let k = 0; k < centroids.length; k++){
            centroids[k].data[0] = 0;
            centroids[k].data[1] = 0;
            centroids[k].data[2] = 0;
            for (let i = 0; i < centroid_elements[k].length; i++){
                centroids[k].data[0] += centroid_elements[k][i].data[0];
                centroids[k].data[1] += centroid_elements[k][i].data[1];
                centroids[k].data[2] += centroid_elements[k][i].data[2];
            }
            centroids[k].data[0] = centroids[k].data[0]/centroid_elements[k].length;
            centroids[k].data[1] = centroids[k].data[1]/centroid_elements[k].length;
            centroids[k].data[2] = centroids[k].data[2]/centroid_elements[k].length;
        }
        
    }

    // copy image over and add color palette
    const newImg = image.createImg(image.width + 50, image.height);

    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++){
            newImg.setPixel(x,y, image.getPixel(x,y));
        }
    } 

    // Adding a palette of size proportionate to the number of pixels.
    const percentage = [];
    let num_pixels = image.height*image.width;
    for(let k = 0; k < centroid_elements.length; k++){
        percentage[k] = centroid_elements[k].length/num_pixels;
    }

    let ybranch = 0;
    let current_color = 0;
    for (let y = 0; y < image.height; y++){
        for (let x = image.width; x < newImg.width; x++) {
        
            if((y- ybranch)/image.height > percentage[current_color]){
                ybranch = y;
                current_color += 1;
                newImg.setPixel(x,y, centroids[current_color]);
            }
            else{
                newImg.setPixel(x,y, centroids[current_color]);
            }
        }
    } 

    // ----------- STUDENT CODE END ------------
    //Gui.alertOnce ('paletteFilter is not implemented yet');
    return newImg;
};

// Read the following paper and implement your own "painter":
//      http://mrl.nyu.edu/publications/painterly98/hertzmann-siggraph98.pdf
Filters.paintFilter = function(image, value) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 59 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('paintFilter is not implemented yet');
    return image;
};

/*
* Read this paper for background on eXtended Difference-of-Gaussians:
*      http://www.cs.princeton.edu/courses/archive/spring19/cos426/papers/Winnemoeller12.pdf
* Read this paper for an approach that develops a flow field based on a bilateral filter
*      http://www.cs.princeton.edu/courses/archive/spring19/cos426/papers/Kang09.pdf
*/
Filters.xDoGFilter = function(image, value) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 70 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('xDoGFilter is not implemented yet');
    return image;
};

// You can use this filter to do whatever you want, for example
// trying out some new idea or implementing something for the
// art contest.
// Currently the 'value' argument will be 1 or whatever else you set
// it to in the URL. You could use this value to switch between
// a bunch of different versions of your code if you want to
// code up a bunch of different things for the art contest.
Filters.customFilter = function(image, value) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 0 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('customFilter is not implemented yet');
    return image;
};
