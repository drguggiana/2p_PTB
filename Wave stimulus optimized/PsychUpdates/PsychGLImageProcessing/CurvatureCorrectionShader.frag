/* Shader for application of a curvature correction of textures during drawing.
 *
 * (C)2012 Pieter Goltstein (Inspired by Labrigger; Spencer Smith)
 *
 */

#extension GL_ARB_texture_rectangle : enable

uniform int           doFilter;
uniform sampler2DRect Image;
uniform vec2          inSize;
uniform vec2          inOffset;
uniform vec2          outSize;

/* Width of the projection area at a distance R in cm */
uniform vec2          ScrSz;

/* Radius of the used spherical screen in cm */
uniform float         ScrDist;

void main()
{
    vec4 texcolor, tl, tr, bl, br;

    /* Get wanted output pixel coordinates for which we should remap:  */
    /* texoutpos is the location of the output pixel for which we need */
    /* to compute the final remapped color value                       */
    vec2 texoutpos = (gl_TexCoord[0].st) - vec2(0.5, 0.5);

    /* Apply 2D offset to output pixel location, as well as a 2D rescaling: */
    /* This to account for situations where not the full display area of    */
    /* display device is used. */
    texoutpos = clamp(texoutpos, vec2(0.0), outSize) / outSize;

    float ScreenDistance = ScrDist;
    float ScreenHalfWidthCm = ScrSz.x/2;
    float ScreenHalfHeightCm = ScrSz.y/2;
    float ScreenHalfWidthRad = atan( ScreenHalfWidthCm, ScreenDistance );
    float ScreenHalfHeightRad = atan( ScreenHalfHeightCm, ScreenDistance );

    texoutpos.x = ((texoutpos.x * 2.0) - 1.0) * ScreenHalfWidthCm;
    texoutpos.y = ((texoutpos.y * 2.0) - 1.0) * ScreenHalfHeightCm;

    float Diam = sqrt(pow(abs(texoutpos.x), 2.0) + pow(abs(ScreenDistance), 2.0));
    texoutpos.x = (atan( texoutpos.x, ScreenDistance )) / ScreenHalfWidthRad;
    texoutpos.y = (atan( texoutpos.y, Diam )) / ScreenHalfHeightRad;

    texoutpos = (texoutpos+1.0)*0.5;

    /* Apply some rescaling to a total input image size of inSize, in case */
    /* that not the full size of the input image is used:                  */
    vec2 texinpos = inOffset + (texoutpos * inSize);

    /* texinpos now contains the location of the final input image pixel to read   */
    /* for creation of the final output color pixel. This value may be fractional. */

    /* Need to filter input image texture ourselves? */
    if (doFilter > 0) {
        /* Need to do our own bilinear filtering: Get colors of four nearest */
        /* neighbour texels for given non-integral, fractional 'texinpos':   */
        tl=texture2DRect(Image, floor(texinpos));
        tr=texture2DRect(Image, floor(texinpos) + vec2(1.0, 0.0));
        bl=texture2DRect(Image, floor(texinpos) + vec2(0.0, 1.0));
        br=texture2DRect(Image, floor(texinpos) + vec2(1.0, 1.0));

        /* Perform weighted linear interpolation -- bilinear interpolation of the 4: */
        tl=mix(tl,tr,fract(texinpos.x));
        bl=mix(bl,br,fract(texinpos.x));
        texcolor = mix(tl, bl, fract(texinpos.y));
    }
    else
    {
        /* Standard case: Hardware does bilinear filtering: */
        texcolor = texture2DRect(Image, texinpos);
    }

    /* Set all undefined (outside input image) pixels to black background color. */
    if (texinpos.x < 0.0 || texinpos.x > inSize.x || texinpos.y < 0.0 || texinpos.y > inSize.y) {
        gl_FragColor = vec4(0.0);
    }
    else
    {
        /* Multiply filtered texcolor with incoming fragment color (GL_MODULATE emulation): */
        /* Assign result as output fragment color: */
        gl_FragColor = texcolor * gl_Color;
    }
}
