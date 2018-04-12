/*------------------------------------------------------------------------
 *  Copyright 2007-2010 (c) Jeff Brown <spadix@users.sourceforge.net>
 *
 *  This file is part of the ZBar Bar Code Reader.
 *
 *  The ZBar Bar Code Reader is free software; you can redistribute it
 *  and/or modify it under the terms of the GNU Lesser Public License as
 *  published by the Free Software Foundation; either version 2.1 of
 *  the License, or (at your option) any later version.
 *
 *  The ZBar Bar Code Reader is distributed in the hope that it will be
 *  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser Public License
 *  along with the ZBar Bar Code Reader; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 *  Boston, MA  02110-1301  USA
 *
 *  http://sourceforge.net/projects/zbar
 *------------------------------------------------------------------------*/

#include "error.h"
#include "image.h"
#include "refcnt.h"

zbar_image_t *zbar_image_create ()
{
    zbar_image_t *img = calloc(1, sizeof(zbar_image_t));
    _zbar_refcnt_init();
    _zbar_image_refcnt(img, 1);
    img->srcidx = -1;
    return(img);
}

void _zbar_image_free (zbar_image_t *img)
{
    if(img->syms) {
        zbar_symbol_set_ref(img->syms, -1);
        img->syms = NULL;
    }
    free(img);
}

void zbar_image_destroy (zbar_image_t *img)
{
    _zbar_image_refcnt(img, -1);
}

void zbar_image_ref (zbar_image_t *img,
                     int refs)
{
    _zbar_image_refcnt(img, refs);
}

unsigned long zbar_image_get_format (const zbar_image_t *img)
{
    return(img->format);
}

unsigned zbar_image_get_sequence (const zbar_image_t *img)
{
    return(img->seq);
}

unsigned zbar_image_get_width (const zbar_image_t *img)
{
    return(img->width);
}

unsigned zbar_image_get_height (const zbar_image_t *img)
{
    return(img->height);
}

void zbar_image_get_size (const zbar_image_t *img,
                          unsigned *w,
                          unsigned *h)
{
    if(w) *w = img->width;
    if(h) *h = img->height;
}

void zbar_image_get_crop (const zbar_image_t *img,
                          unsigned *x,
                          unsigned *y,
                          unsigned *w,
                          unsigned *h)
{
    if(x) *x = img->crop_x;
    if(y) *y = img->crop_y;
    if(w) *w = img->crop_w;
    if(h) *h = img->crop_h;
}

const void *zbar_image_get_data (const zbar_image_t *img)
{
    return(img->data);
}

unsigned long zbar_image_get_data_length (const zbar_image_t *img)
{
    return(img->datalen);
}

void zbar_image_set_format (zbar_image_t *img,
                            unsigned long fmt)
{
    img->format = fmt;
}

void zbar_image_set_sequence (zbar_image_t *img,
                              unsigned seq)
{
    img->seq = seq;
}

void zbar_image_set_size (zbar_image_t *img,
                          unsigned w,
                          unsigned h)
{
    img->crop_x = img->crop_y = 0;
    img->width = img->crop_w = w;
    img->height = img->crop_h = h;
}

void zbar_image_set_crop (zbar_image_t *img,
                          unsigned x,
                          unsigned y,
                          unsigned w,
                          unsigned h)
{
    unsigned img_w = img->width;
    if(x > img_w) x = img_w;
    if(x + w > img_w) w = img_w - x;
    img->crop_x = x;
    img->crop_w = w;

    unsigned img_h = img->height;
    if(y > img_h) y = img_h;
    if(y + h > img_h) h = img_h - y;
    img->crop_y = y;
    img->crop_h = h;
}

inline void zbar_image_free_data (zbar_image_t *img)
{
    if(!img)
        return;
    if(img->src) {
        zbar_image_t *newimg;
        /* replace video image w/new copy */
        assert(img->refcnt); /* FIXME needs lock */
        newimg = zbar_image_create();
        memcpy(newimg, img, sizeof(zbar_image_t));
        /* recycle video image */
        newimg->cleanup(newimg);
        /* detach old image from src */
        img->cleanup = NULL;
        img->src = NULL;
        img->srcidx = -1;
    }
    else if(img->cleanup && img->data) {
        if(img->cleanup != zbar_image_free_data) {
            /* using function address to detect this case is a bad idea;
             * windows link libraries add an extra layer of indirection...
             * this works around that problem (bug #2796277)
             */
            zbar_image_cleanup_handler_t *cleanup = img->cleanup;
            img->cleanup = zbar_image_free_data;
            cleanup(img);
        }
        else
            free((void*)img->data);
    }
    img->data = NULL;
}

void zbar_image_set_data (zbar_image_t *img,
                          const void *data,
                          unsigned long len,
                          zbar_image_cleanup_handler_t *cleanup)
{
    zbar_image_free_data(img);
    img->data = data;
    img->datalen = len;
    img->cleanup = cleanup;
}

void zbar_image_set_userdata (zbar_image_t *img,
                              void *userdata)
{
    img->userdata = userdata;
}

void *zbar_image_get_userdata (const zbar_image_t *img)
{
    return(img->userdata);
}

zbar_image_t *zbar_image_copy (const zbar_image_t *src)
{
    zbar_image_t *dst = zbar_image_create();
    dst->format = src->format;
    _zbar_image_copy_size(dst, src);
    dst->datalen = src->datalen;
    dst->data = malloc(src->datalen);
    assert(dst->data);
    memcpy((void*)dst->data, src->data, src->datalen);
    dst->cleanup = zbar_image_free_data;
    return(dst);
}

zbar_image_t *zbar_image_empty_copy (const zbar_image_t *src)
{
    zbar_image_t *dst = zbar_image_create();
    dst->format = src->format;
    _zbar_image_copy_size(dst, src);
    dst->datalen = src->datalen;
    dst->data = malloc(src->datalen);
    assert(dst->data);
    dst->cleanup = zbar_image_free_data;
    return(dst);
}


zbar_image_t *zbar_image_diff (const zbar_image_t *src, const zbar_image_t *sub, zbar_image_t *dst)
{
    const uint8_t *scan = src->data;
    const uint8_t *scan2 = sub->data;
    if (!dst) {
        dst = zbar_image_empty_copy(src);
    }
    uint8_t* put = (uint8_t*)dst->data;
    const uint8_t *end = scan + src->width * src->height;
    do {
        int16_t v = *scan++ - *scan2++;
        v &= ~(v >> 16);                                // clamp negative to 0
        v |= (0xff - v) >> 16;                          // if >= 0x100, set all bits
        *put++ = (uint8_t)v;
    } while (scan != end);
    return dst;
}

// adjust brightness (gain) and contrast (bias)
zbar_image_t *zbar_image_adjust(const zbar_image_t *src, zbar_image_t *dst, double gain, double bias)
{
    const uint8_t *scan = src->data;
    if (!dst) {
        dst = zbar_image_empty_copy(src);
    }
    uint8_t* put = (uint8_t*)dst->data;
    const uint8_t *end = scan + src->width * src->height;
    int gf = (int)(gain * 0x100);
    int bf = (int)(bias * 0x10000);
    do {
        int16_t v = *scan++;
        v = ((v - 0x80) * gf + 0x80 + bf) >> 8;
        v &= ~(v >> 16);                                // clamp negative to 0
        v |= (0xff - v) >> 16;                          // if >= 0x100, set all bits
        *put++ = (uint8_t)v;
    } while (scan != end);
    return dst;
}


static void zbar_image_filter33 (const zbar_image_t *src, zbar_image_t *dst, const int16_t* kernel)
{
    uint8_t* put = (uint8_t*)dst->data;
    unsigned w = src->width;
    const uint8_t *scan = src->data;
    uint8_t *end = put + src->width * (src->height - 1) - 1;
    memcpy(put, scan, w + 1);
    put += src->width + 1;
    scan += src->width + 1;
    do {
        int16_t v = (
            scan[-w-1] * (int)kernel[1] +
            scan[-w  ] * (int)kernel[2] +
            scan[-w+1] * (int)kernel[3] +
            scan[  -1] * (int)kernel[4] +
            scan[   0] * (int)kernel[5] +
            scan[   1] * (int)kernel[6] +
            scan[ w-1] * (int)kernel[7] +
            scan[ w  ] * (int)kernel[8] +
            scan[ w+1] * (int)kernel[9]) >> 8;
        v &= ~(v >> 16);                                // clamp negative to 0
        v |= (0xff - v) >> 16;                          // if >= 0x100, set all bits
        *put++ = (uint8_t)v;
        ++scan;
    } while (put != end);
    memcpy(put, scan, w + 1);
}

#include <math.h>
static int *boxesForGauss(double sigma, int n, int *buffer)  // standard deviation, number of boxes
{
    double wIdeal = sqrt((12*sigma*sigma/n)+1);  // Ideal averaging filter width
    int wl = floor(wIdeal);
    wl -= (wl & 1) ^ 1;
    int wu = wl + 2;

    double mIdeal = (12*sigma*sigma - n*wl*wl - 4*n*wl - 3*n)/(-4*wl - 4);
    double m = round(mIdeal);
    // var sigmaActual = Math.sqrt( (m*wl*wl + (n-m)*wu*wu - n)/12 );

    for (int i = 0; i < n; i++) {
        buffer[i] = i < m ? wl : wu;
    }
    return buffer;
}

static void boxBlurT4(const uint8_t* src, uint8_t* dst, int w, int h, int r)
{
    int iarr = (1.0 / (r + r + 1)) * 0x10000;
    for (int i = 0; i < w; i++) {
        int ti = i, li = ti, ri = ti + r * w;
        int fv = src[ti], lv = src[ti + w * (h - 1)], val = (r + 1) * fv;
        const uint8_t *scan = src + ti;
        const uint8_t *end = scan + r * w;
        do {
            val += *scan;
            scan += w;
        } while (scan != end);
        scan = src + ri;
        end = scan + (r + 1) * w;
        uint8_t *put = dst + ti;
        do {
            val += *scan - fv;
            *put = (val * iarr + 0x8000) >> 16;
            scan += w;
            put += w;
        } while (scan != end);
        end = src + ti + (h - r) * w;
        const uint8_t *prev = src + li;
        do {
            val += *scan - *prev;
            *put = (val * iarr + 0x8000) >> 16;
            scan += w;
            prev += w;
            put += w;
        } while (scan != end);
        end = dst + ti + (h * w);
        do {
            val += lv - *prev;
            *put = (val * iarr + 0x8000) >> 16;
            prev += w;
            put += w;
        } while (put != end);
    }
}

static void boxBlurH4(const uint8_t* src, uint8_t* dst, int w, int h, int r)
{
    int iarr = (1.0 / (r + r + 1)) * 0x10000;
    for (int i = 0; i < h; i++) {
        int ti = i * w, li = ti, ri = ti + r;
        int fv = src[ti], lv = src[ti + w - 1], val = (r + 1) * fv;
        const uint8_t *scan = src + ti;
        const uint8_t *end = scan + r;
        do {
            val += *scan++;
        } while (scan != end);
        scan = src + ri;
        end = scan + r + 1;
        uint8_t *put = dst + ti;
        do {
            val += *scan++ - fv;
            *put++ = (val * iarr + 0x8000) >> 16;
        } while (scan != end);
        end = src + ti + w - r;
        const uint8_t *prev = src + li;
        do {
            val += *scan++ - *prev++;
            *put++ = (val * iarr + 0x8000) >> 16;
        } while (scan != end);
        end = src + ti + w;
        do {
            val += lv - *prev++;
            *put++ = (val * iarr + 0x8000) >> 16;
        } while (++scan != end);
    }
}

static void boxBlur4(const uint8_t* src, uint8_t* dst, int w, int h, int r)
{
    memcpy(dst, src, w * h);
    boxBlurH4(src, dst, w, h, r);
    boxBlurT4(src, dst, w, h, r);
}

zbar_image_t *zbar_image_blur (const zbar_image_t *src, zbar_image_t *dst, int r)
{
    int bxs[3];
    if (!dst) {
        dst = zbar_image_empty_copy(src);
    }
    boxesForGauss(r, 3, bxs);
    boxBlur4(src->data, (uint8_t*)dst->data, src->width, src->height, (bxs[0] - 1)/2);
    boxBlur4(dst->data, (uint8_t*)dst->data, src->width, src->height, (bxs[1] - 1)/2);
    boxBlur4(dst->data, (uint8_t*)dst->data, src->width, src->height, (bxs[2] - 1)/2);
    return dst;
}

zbar_image_t *zbar_image_filter (const zbar_image_t *src, zbar_image_t *dst, const int16_t *kernel)
{
    if (!dst) {
        dst = zbar_image_empty_copy(src);
    }
    if (kernel[0] == 9) {
        zbar_image_filter33(src, dst, kernel);
    } else if (kernel[0] == 25) {
//        zbar_image_filter55(src, dst, kernel);
    }
    return(dst);
}

zbar_image_t *zbar_image_unsharp_mask (const zbar_image_t *src, double amount, double radius, double threshold)
{
    zbar_image_t* blur;
    zbar_image_t* mask;
    zbar_image_t* contrast;
    /*
    int16_t contrastFilter[10]={
        9,
        -0x0100, -0x0100, -0x0100,
        -0x0100,  0x0800, -0x0100,
        -0x0100, -0x0100, -0x0100
    };
     */
    blur = zbar_image_blur(src, NULL, (int)(radius + 0.99));
    //contrast = zbar_image_filter(src, NULL, contrastFilter);
    contrast = zbar_image_adjust(src, NULL, amount, 0);
    mask = zbar_image_diff(src, blur, blur);
    const uint8_t *scan = src->data;
    const uint8_t *end = src->data + src->width * src->height;
    const uint8_t *ctscan = contrast->data;
    uint8_t *uscan = (uint8_t *)mask->data;
    int tf = threshold * 0x100;
    do {
        int16_t v;
        v = *scan++;
        int d = v - *ctscan++;
        int pct = *uscan;
        int delta = (d * pct) >> 8;
        if (abs(delta) > tf) {
            v += delta;
            v &= ~(v >> 16);                                // clamp negative to 0
            v |= (0xff - v) >> 16;                          // if >= 0x100, set all bits
        }
        *uscan = (uint8_t)v;
        ++uscan;
    } while (scan != end);
    zbar_image_destroy(contrast);
    return mask;
}

const zbar_symbol_set_t *zbar_image_get_symbols (const zbar_image_t *img)
{
    return(img->syms);
}

void zbar_image_set_symbols (zbar_image_t *img,
                             const zbar_symbol_set_t *syms)
{
    if(syms)
        zbar_symbol_set_ref(syms, 1);
    if(img->syms)
        zbar_symbol_set_ref(img->syms, -1);
    img->syms = (zbar_symbol_set_t*)syms;
}

const zbar_symbol_t *zbar_image_first_symbol (const zbar_image_t *img)
{
    return((img->syms) ? img->syms->head : NULL);
}

typedef struct zimg_hdr_s {
    uint32_t magic, format;
    uint16_t width, height;
    uint32_t size;
} zimg_hdr_t;

int zbar_image_write (const zbar_image_t *img,
                      const char *filebase)
{
    int len = strlen(filebase) + 16;
    char *filename = malloc(len);
    int n = 0, rc = 0;
    FILE *f;
    zimg_hdr_t hdr;
    strcpy(filename, filebase);
    if((img->format & 0xff) >= ' ')
        n = snprintf(filename, len, "%s.%.4s.zimg",
                     filebase, (char*)&img->format);
    else
        n = snprintf(filename, len, "%s.%08" PRIx32 ".zimg",
                     filebase, img->format);
    assert(n < len - 1);
    filename[len - 1] = '\0';

    zprintf(1, "dumping %.4s(%08" PRIx32 ") image to %s\n",
            (char*)&img->format, img->format, filename);

    f = fopen(filename, "w");
    if(!f) {
#ifdef HAVE_ERRNO_H
        rc = errno;
        zprintf(1, "ERROR opening %s: %s\n", filename, strerror(rc));
#else
        rc = 1;
#endif
        goto error;
    }

    hdr.magic = 0x676d697a;
    hdr.format = img->format;
    hdr.width = img->width;
    hdr.height = img->height;
    hdr.size = img->datalen;

    if(fwrite(&hdr, sizeof(hdr), 1, f) != 1 ||
       fwrite(img->data, 1, img->datalen, f) != img->datalen) {
#ifdef HAVE_ERRNO_H
        rc = errno;
        zprintf(1, "ERROR writing %s: %s\n", filename, strerror(rc));
#else
        rc = 1;
#endif
        fclose(f);
        goto error;
    }

    rc = fclose(f);

error:
    free(filename);
    return(rc);
}

#ifdef DEBUG_SVG
# include <png.h>

int zbar_image_write_png (const zbar_image_t *img,
                          const char *filename)
{
    int rc = -1;
    FILE *file = NULL;
    png_struct *png = NULL;
    png_info *info = NULL;
    const uint8_t **rows = NULL;

    rows = malloc(img->height * sizeof(*rows));
    if(!rows)
        goto done;

    rows[0] = img->data;
    int y;
    for(y = 1; y < img->height; y++)
        rows[y] = rows[y - 1] + img->width;

    file = fopen(filename, "wb");
    if(!file)
        goto done;

    png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if(!png)
        goto done;

    info = png_create_info_struct(png);
    if(!info)
        goto done;

    if(setjmp(png_jmpbuf(png)))
        goto done;

    png_init_io(png, file);
    png_set_compression_level(png, 9);
    png_set_IHDR(png, info, img->width, img->height, 8, PNG_COLOR_TYPE_GRAY,
                 PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
                 PNG_FILTER_TYPE_DEFAULT);

    png_set_rows(png, info, (void*)rows);
    png_write_png(png, info, PNG_TRANSFORM_IDENTITY, NULL);

    png_write_end(png,info);
    rc = 0;

done:
    if(png)
        png_destroy_write_struct(&png, &info);
    if(rows)
        free(rows);
    if(file)
        fclose(file);
    return(rc);
}

#endif
