/*! \file local_align.c
 *  \brief Calculates the local sequence alignment by Smith-Waterman algorithm
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: September 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <emmintrin.h>
#include "ddradseq.h"

#define KSW_XBYTE  0x10000
#define KSW_XSTOP  0x20000
#define KSW_XSUBO  0x40000
#define KSW_XSTART 0x80000
#define MINUS_INF -0x40000000

const ALIGN_RESULT g_defr = { 0, -1, -1, -1, -1, -1, -1 };

ALIGN_QUERY *
align_init (int size, int qlen, const unsigned char *query, int m,
            const char *mat)
{
    int slen = 0;
    int a = 0;
    int tmp = 0;
    int p = 0;
    ALIGN_QUERY *q = NULL;

    size = (size > 1) ? 2 : 1;

    /* Number of values per __m128i */
    p = 8 * (3 - size);

    /* Segmented length */
    slen = (qlen + p - 1) / p;

    /* Allocate memory for query profile */
    if ((q = malloc (sizeof (ALIGN_QUERY) + 256 + 16 * slen * (m + 4))) == NULL)
        {
            fprintf (stderr, "[libfasta:%s:%d] Error: cannot allocate "
                     "memory for query profile.\n", __func__, __LINE__);
            return NULL;
        }

    /* Align memory */
    q->qp = (__m128i*)(((size_t)q + sizeof (ALIGN_QUERY) + 15) >> 4 << 4);
    q->H0 = q->qp + slen * m;
    q->H1 = q->H0 + slen;
    q->E  = q->H1 + slen;
    q->Hmax = q->E + slen;
    q->slen = slen;
    q->qlen = qlen;
    q->size = size;

    /* Compute shift */
    tmp = m * m;

    /* Find the minimum and maximum score */
    for (a = 0, q->shift = 127, q->mdiff = 0; a < tmp; a++)
        {
            if (mat[a] < (char)q->shift)
                {
                    q->shift = mat[a];
                }

            if (mat[a] > (char)q->mdiff)
                {
                    q->mdiff = mat[a];
                }
        }

    q->max = q->mdiff;

    /* NB: q->shift is uint8_t */
    q->shift = 256 - q->shift;

    /* Difference between the min and max scores */
    q->mdiff += q->shift;

    /* An example: p=8, qlen=19, slen=3 and segmentation:
     * {{0,3,6,9,12,15,18,-1},{1,4,7,10,13,16,-1,-1},{2,5,8,11,14,17,-1,-1}} */
    if (size == 1)
        {
            char *t = (char*)q->qp;

            for (a = 0; a < m; a++)
                {
                    int i = 0;
                    int k = 0;
                    int nlen = slen * p;
                    const char *ma = mat + a * m;

                    for (i = 0; i < slen; ++i)
                        {
                            /* p iterations */
                            for (k = i; k < nlen; k += slen)
                                {
                                    *t++ = ((k >= qlen) ? 0 : ma[query[k]]) + q->shift;
                                }
                        }
                }
        }
    else
        {
            short int *t = (short int*)q->qp;

            for (a = 0; a < m; a++)
                {
                    int i = 0;
                    int k = 0;
                    int nlen = slen * p;
                    const char *ma = mat + a * m;

                    for (i = 0; i < slen; i++)
                        {
                            /* p iterations */
                            for (k = i; k < nlen; k += slen)
                                {
                                    *t++ = ((k >= qlen) ? 0 : ma[query[k]]);
                                }
                        }
                }
        }

    return q;
}

ALIGN_RESULT
ksw_u8 (ALIGN_QUERY *q, int tlen, const unsigned char *target, int _gapo, int _gape, int xtra)
{
    int slen = 0;
    int i = 0;
    int m_b = 0;
    int n_b = 0;
    int te = -1;
    int gmax = 0;
    int minsc = 0;
    int endsc = 0;
    uint64_t *b;
    __m128i zero;
    __m128i gapoe;
    __m128i gape;
    __m128i shift;
    __m128i *H0;
    __m128i *H1;
    __m128i *E;
    __m128i *Hmax;
    ALIGN_RESULT r;

#define __max_16(ret, xx) do { \
        (xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 8)); \
        (xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 4)); \
        (xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 2)); \
        (xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 1)); \
        (ret) = _mm_extract_epi16((xx), 0) & 0x00ff; \
    } while (0)

    /* Initialization */
    r = g_defr;
    minsc = (xtra & KSW_XSUBO) ? xtra & 0xffff : 0x10000;
    endsc = (xtra & KSW_XSTOP) ? xtra & 0xffff : 0x10000;
    m_b = n_b = 0;
    b = 0;
    zero = _mm_set1_epi32 (0);
    gapoe = _mm_set1_epi8 (_gapo + _gape);
    gape = _mm_set1_epi8 (_gape);
    shift = _mm_set1_epi8 (q->shift);
    H0 = q->H0;
    H1 = q->H1;
    E = q->E;
    Hmax = q->Hmax;
    slen = q->slen;

    for (i = 0; i < slen; i++)
        {
            _mm_store_si128 (E + i, zero);
            _mm_store_si128 (H0 + i, zero);
            _mm_store_si128 (Hmax + i, zero);
        }

    /* Core loop */
    for (i = 0; i < tlen; i++)
        {
            int j = 0;
            int k = 0;
            int cmp = 0;
            int imax = 0;
            __m128i e;
            __m128i h;
            __m128i f = zero;
            __m128i max = zero;
            __m128i *S = q->qp + target[i] * slen;

            /* h={2,5,8,11,14,17,-1,-1} in the above example */
            h = _mm_load_si128 (H0 + slen - 1);

            /* h=H(i-1,-1); << instead of >> because x64 is little-endian */
            h = _mm_slli_si128 (h, 1);

            for (j = 0; LIKELY(j < slen); j++)
                {
                    /* SW cells are computed in the following order:
                     *   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
                     *   E(i+1,j) = max{H(i,j)-q, E(i,j)-r}
                     *   F(i,j+1) = max{H(i,j)-q, F(i,j)-r}
                     */

                    /* Compute H'(i,j); note that at the beginning, h=H'(i-1,j-1) */
                    h = _mm_adds_epu8 (h, _mm_load_si128(S + j));

                    /* h=H'(i-1,j-1)+S(i,j) */
                    h = _mm_subs_epu8 (h, shift);

                    /* e=E'(i,j) */
                    e = _mm_load_si128 (E + j);
                    h = _mm_max_epu8 (h, e);

                    /* h=H'(i,j) */
                    h = _mm_max_epu8 (h, f);

                    /* Set max */
                    max = _mm_max_epu8 (max, h);

                    /* Save to H'(i,j) */
                    _mm_store_si128 (H1 + j, h);

                    /* Now compute E'(i+1,j) */
                    /* h=H'(i,j)-gapo */
                    h = _mm_subs_epu8 (h, gapoe);

                    /* e=E'(i,j)-gape */
                    e = _mm_subs_epu8 (e, gape);

                    /* e=E'(i+1,j) */
                    e = _mm_max_epu8 (e, h);

                    /* Save to E'(i+1,j) */
                    _mm_store_si128 (E + j, e);

                    /* Now compute F'(i,j+1) */
                    f = _mm_subs_epu8 (f, gape);
                    f = _mm_max_epu8 (f, h);

                    /* get H'(i-1,j) and prepare for the next j */
                    /* h=H'(i-1,j) */
                    h = _mm_load_si128 (H0 + j);
                }

            /* NB: we do not need to set E(i,j) as we disallow */
            /* adjecent insertion and then deletion */
            /* this block mimics SWPS3; NB: H(i,j) updated in the */
            /* lazy-F loop cannot exceed max */
            for (k = 0; LIKELY(k < 16); k++)
                {
                    f = _mm_slli_si128 (f, 1);

                    for (j = 0; LIKELY(j < slen); j++)
                        {
                            h = _mm_load_si128 (H1 + j);

                            /* h=H'(i,j) */
                            h = _mm_max_epu8 (h, f);
                            _mm_store_si128 (H1 + j, h);
                            h = _mm_subs_epu8 (h, gapoe);
                            f = _mm_subs_epu8 (f, gape);
                            cmp = _mm_movemask_epi8 (_mm_cmpeq_epi8 (_mm_subs_epu8 (f, h), zero));

                            if (UNLIKELY(cmp == 0xffff))
                                {
                                    goto end_loop16;
                                }
                        }
                }
end_loop16:
            /* imax is the maximum number in max */
            __max_16(imax, max);

            /* Write the b array; this condition adds branching unfornately */
            if (imax >= minsc)
                {
                    /* Then append */
                    if ((n_b == 0) || ((int)b[n_b-1] + 1 != i))
                        {
                            if (n_b == m_b)
                                {
                                    m_b = m_b ? m_b << 1 : 8;
                                    if ((b = realloc (b, 8 * m_b)) == NULL)
                                        {
                                            fprintf (stderr, "[libfasta:%s:%d] Error: cannot allocate "
                                                     "memory for b.\n", __func__, __LINE__);
                                            return r;
                                        }
                                }
                            b[n_b++] = (uint64_t)imax << 32 | i;
                        }
                    else if ((int)(b[n_b-1] >> 32) < imax)
                        {
                            /* Modify the last */
                            b[n_b-1] = (uint64_t)imax << 32 | i;
                        }
                }

            if (imax > gmax)
                {
                    gmax = imax;

                    /* te is the end position on the target */
                    te = i;

                    /* Keep the H1 vector */
                    for (j = 0; LIKELY(j < slen); j++)
                        {
                            _mm_store_si128 (Hmax + j, _mm_load_si128 (H1 + j));
                        }

                    if ((gmax + q->shift >= 255) || (gmax >= endsc))
                        {
                            break;
                        }
                }
            S = H1;
            H1 = H0;

            /* Swap H0 and H1 */
            H0 = S;
        }

    r.score = (gmax + q->shift < 255) ? gmax : 255;
    r.target_end = te;

    /* Get a->query_end, the end of query match */
    /* find the 2nd best score */
    if (r.score != 255)
        {
            int max = -1;
            int low = 0;
            int high = 0;
            int qlen = slen * 16;
            unsigned char *t = (unsigned char *)Hmax;

            for (i = 0; i < qlen; i++, t++)
                {
                    if ((int)*t > max)
                        {
                            max = *t;
                            r.query_end = i / 16 + i % 16 * slen;
                        }
                }

            if (b)
                {
                    i = (r.score + q->max - 1) / q->max;
                    low = te - i;
                    high = te + i;

                    for (i = 0; i < n_b; i++)
                        {
                            int e = (int)b[i];

                            if ((e < low || e > high) && (int)(b[i] >> 32) > r.score2)
                                {
                                    r.score2 = b[i] >> 32;
                                    r.target_end2 = e;
                                }
                        }
                }
        }

    free (b);

    return r;
}

ALIGN_RESULT
ksw_i16 (ALIGN_QUERY *q, int tlen, const unsigned char *target,
         int _gapo, int _gape, int xtra)
{
    int slen = 0;
    int i = 0;
    int m_b = 0;
    int n_b = 0;
    int te = -1;
    int gmax = 0;
    int minsc = 0;
    int endsc = 0;
    uint64_t *b = 0;
    __m128i zero;
    __m128i gapoe;
    __m128i gape;
    __m128i *H0;
    __m128i *H1;
    __m128i *E;
    __m128i *Hmax;
    ALIGN_RESULT r;

#define __max_8(ret, xx) do { \
        (xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 8)); \
        (xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 4)); \
        (xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 2)); \
        (ret) = _mm_extract_epi16((xx), 0); \
    } while (0)

    /* Initialization */
    r = g_defr;
    minsc = (xtra & KSW_XSUBO) ? xtra & 0xffff : 0x10000;
    endsc = (xtra & KSW_XSTOP) ? xtra & 0xffff : 0x10000;
    m_b = n_b = 0;
    b = 0;
    zero = _mm_set1_epi32 (0);
    gapoe = _mm_set1_epi16 (_gapo + _gape);
    gape = _mm_set1_epi16 (_gape);
    H0 = q->H0;
    H1 = q->H1;
    E = q->E;
    Hmax = q->Hmax;
    slen = q->slen;

    for (i = 0; i < slen; ++i)
        {
            _mm_store_si128 (E + i, zero);
            _mm_store_si128 (H0 + i, zero);
            _mm_store_si128 (Hmax + i, zero);
        }

    /* Core loop */
    for (i = 0; i < tlen; ++i)
        {
            int j = 0;
            int k = 0;
            int imax = 0;
            __m128i e;
            __m128i h;
            __m128i f = zero;
            __m128i max = zero;
            __m128i *S = q->qp + target[i] * slen;

            /* h={2,5,8,11,14,17,-1,-1} in the above example */
            h = _mm_load_si128 (H0 + slen - 1);
            h = _mm_slli_si128 (h, 2);

            for (j = 0; LIKELY(j < slen); j++)
                {
                    h = _mm_adds_epi16 (h, *S++);
                    e = _mm_load_si128 (E + j);
                    h = _mm_max_epi16 (h, e);
                    h = _mm_max_epi16 (h, f);
                    max = _mm_max_epi16 (max, h);
                    _mm_store_si128 (H1 + j, h);
                    h = _mm_subs_epu16 (h, gapoe);
                    e = _mm_subs_epu16 (e, gape);
                    e = _mm_max_epi16 (e, h);
                    _mm_store_si128 (E + j, e);
                    f = _mm_subs_epu16 (f, gape);
                    f = _mm_max_epi16 (f, h);
                    h = _mm_load_si128 (H0 + j);
                }

            for (k = 0; LIKELY(k < 16); k++)
                {
                    f = _mm_slli_si128 (f, 2);
                    for (j = 0; LIKELY(j < slen); j++)
                        {
                            h = _mm_load_si128 (H1 + j);
                            h = _mm_max_epi16 (h, f);
                            _mm_store_si128 (H1 + j, h);
                            h = _mm_subs_epu16 (h, gapoe);
                            f = _mm_subs_epu16 (f, gape);

                            if (UNLIKELY(!_mm_movemask_epi8 (_mm_cmpgt_epi16 (f, h))))
                                {
                                    goto end_loop8;
                                }
                        }
                }
end_loop8:
            __max_8(imax, max);

            if (imax >= minsc)
                {
                    if ((n_b == 0) || ((int)b[n_b-1] + 1 != i))
                        {
                            if (n_b == m_b)
                                {
                                    m_b = m_b ? m_b << 1 : 8;
                                    if ((b = realloc (b, 8 * m_b)) == NULL)
                                        {
                                            fprintf (stderr, "[libfasta:%s:%d] Error: cannot allocate "
                                                     "memory for b array.\n", __func__, __LINE__);
                                            return r;
                                        }
                                }
                            b[n_b++] = (uint64_t)imax << 32 | i;
                        }
                    else if ((int)(b[n_b-1] >> 32) < imax)
                        {
                            /* Modify the last */
                            b[n_b-1] = (uint64_t)imax << 32 | i;
                        }
                }

            if (imax > gmax)
                {
                    gmax = imax;
                    te = i;

                    for (j = 0; LIKELY(j < slen); j++)
                        {
                            _mm_store_si128 (Hmax + j, _mm_load_si128 (H1 + j));
                        }

                    if (gmax >= endsc)
                        {
                            break;
                        }
                }
            S = H1;
            H1 = H0;
            H0 = S;
        }

    r.score = gmax;
    r.target_end = te;
    {
        int max = -1;
        int low = 0;
        int high = 0;
        int qlen = slen * 8;
        unsigned short int *t = (unsigned short int *)Hmax;

        for (i = 0, r.query_end = -1; i < qlen; i++, t++)
            {
                if ((int)*t > max)
                    {
                        max = *t;
                        r.query_end = i / 8 + i % 8 * slen;
                    }
            }

        if (b)
            {
                i = (r.score + q->max - 1) / q->max;
                low = te - i;
                high = te + i;

                for (i = 0; i < n_b; i++)
                    {
                        int e = (int)b[i];

                        if (((e < low) || (e > high)) && ((int)(b[i] >> 32) > r.score2))
                            {
                                r.score2 = b[i] >> 32;
                                r.target_end2 = e;
                            }
                    }
            }
    }

    free (b);

    return r;
}

static void
revseq (int l, unsigned char *s)
{
    int i = 0;
    int t = 0;

    for (i = 0; i < l >> 1; ++i)
        {
            t = s[i];
            s[i] = s[l - 1 - i];
            s[l - 1 - i] = t;
        }
}

ALIGN_RESULT
local_align (int qlen, unsigned char *query, int tlen,
             unsigned char *target, int m, const char *mat, int gapo,
             int gape, int xtra, ALIGN_QUERY **qry)
{
    int size = 0;
    ALIGN_QUERY *q;
    ALIGN_RESULT r;
    ALIGN_RESULT rr;
    ALIGN_RESULT (*func)(ALIGN_QUERY *, int, const unsigned char *, int, int, int);

    q = (qry && *qry) ? *qry : align_init ((xtra & KSW_XBYTE) ? 1 : 2, qlen, query, m, mat);

    if (qry && (*qry == 0))
        {
            *qry = q;
        }

    func = (q->size == 2) ? ksw_i16 : ksw_u8;
    size = q->size;
    r = func (q, tlen, target, gapo, gape, xtra);

    if (qry == NULL)
        {
            free (q);
        }

    if (((xtra & KSW_XSTART) == 0) || ((xtra & KSW_XSUBO) && (r.score < (xtra & 0xffff))))
        {
            return r;
        }

    revseq (r.query_end + 1, query);

    /* +1 because qe/te points to the exact end */
    /* not the position after the end */
    revseq (r.target_end + 1, target);
    q = align_init (size, r.query_end + 1, query, m, mat);
    rr = func (q, tlen, target, gapo, gape, KSW_XSTOP | r.score);
    revseq (r.query_end + 1, query);
    revseq (r.target_end + 1, target);
    free (q);

    if (r.score == rr.score)
        {
            r.target_begin = r.target_end - rr.target_end;
            r.query_begin = r.query_end - rr.query_end;
        }

    return r;
}

/********************
 *** SW extension ***
 ********************/

typedef struct
{
    int h;
    int e;
} eh_t;

int
align_extend (int qlen, const unsigned char *query, int tlen,
              const unsigned char *target, int m, const char *mat,
              int gapo, int gape, int w, int h0, int *_qle, int *_tle)
{
    int i = 0;
    int j = 0;
    int k = 0;
    int gapoe = gapo + gape;
    int beg = 0;
    int end = 0;
    int max = 0;
    int max_i = 0;
    int max_j = 0;
    int max_gap = 0;
    char *qp = NULL;
    eh_t *eh;

    if (h0 < 0)
        {
            h0 = 0;
        }

    /* Allocate memory */
    if ((qp = malloc (qlen * m)) == NULL)
        {
            fprintf (stderr, "[libfasta:%s:%d] Error: cannot allocate "
                     "memory for the query profile.\n", __func__, __LINE__);
            return -1;
        }
    if ((eh = calloc (qlen + 1, 8)) == NULL)
        {
            fprintf (stderr, "[libfasta:%s:%d] Error: cannot allocate "
                     "memory for eh.\n", __func__, __LINE__);
            return -1;
        }

    /* Generate the query profile */
    for (k = i = 0; k < m; k++)
        {
            const char *p = &mat[k * m];

            for (j = 0; j < qlen; j++)
                {
                    qp[i++] = p[query[j]];
                }
        }

    /* Fill the first row */
    eh[0].h = h0;
    eh[1].h = (h0 > gapoe) ? h0 - gapoe : 0;

    for (j = 2; (j <= qlen) && (eh[j-1].h > gape); j++)
        {
            eh[j].h = eh[j - 1].h - gape;
        }

    /* Adjust $w if it is too large */
    k = m * m;

    /* Get the max score */
    for (i = 0, max = 0; i < k; i++)
        {
            max = (max > mat[i]) ? max : mat[i];
        }

    max_gap = (int)((double)(qlen * max - gapo) / gape + 1.);
    max_gap = (max_gap > 1) ? max_gap : 1;
    w = (w < max_gap) ? w : max_gap;

    /* DP loop */
    max = h0;
    max_i = max_j = -1;
    beg = 0;
    end = qlen;

    for (i = 0; LIKELY(i < tlen); i++)
        {
            int f = 0;
            int h1 = 0;
            int m = 0;
            int mj = -1;
            char *q = &qp[target[i] * qlen];

            /* Compute the first column */
            h1 = h0 - (gapo + gape * (i + 1));

            if (h1 < 0)
                {
                    h1 = 0;
                }

            /* Apply the band and the constraint (if provided) */
            if (beg < i - w)
                {
                    beg = i - w;
                }

            if (end > i + w + 1)
                {
                    end = i + w + 1;
                }

            if (end > qlen)
                {
                    end = qlen;
                }

            for (j = beg; LIKELY(j < end); j++)
                {
                    /* At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
                     * Similar to SSE2-SW, cells are computed in the following order:
                     *   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
                     *   E(i+1,j) = max{H(i,j)-gapo, E(i,j)} - gape
                     *   F(i,j+1) = max{H(i,j)-gapo, F(i,j)} - gape */
                    eh_t *p = &eh[j];
                    /* Get H(i-1,j-1) and E(i-1,j) */
                    int h = p->h;
                    int e = p->e;

                    /* Set H(i,j-1) for the next row */
                    p->h = h1;
                    h += q[j];
                    h = (h > e) ? h : e;
                    h = (h > f) ? h : f;

                    /* Save H(i,j) to h1 for the next column */
                    h1 = h;
                    mj = (m > h) ? mj : j;

                    /* m is stored at eh[mj+1] */
                    m = (m > h) ? m : h;
                    h -= gapoe;
                    h = (h > 0) ? h : 0;
                    e -= gape;

                    /* Computed E(i+1,j) */
                    e = (e > h) ? e : h;

                    /* Save E(i+1,j) for the next row */
                    p->e = e;
                    f -= gape;

                    /* Computed F(i,j+1) */
                    f = (f > h) ? f : h;
                }

            eh[end].h = h1;
            eh[end].e = 0;

            if (m == 0)
                {
                    break;
                }

            if (m > max)
                {
                    max = m;
                    max_i = i;
                    max_j = mj;
                }

            /* Update beg and end for the next round */
            for (j = mj; j >= beg && eh[j].h; j--);
            beg = j + 1;
            for (j = mj + 2; j <= end && eh[j].h; j++);
            end = j;
            /* beg = 0; end = qlen; */ /* uncomment this line for debugging */
        }

    free (eh);
    free (qp);

    if (_qle)
        {
            *_qle = max_j + 1;
        }

    if (_tle)
        {
            *_tle = max_i + 1;
        }

    return max;
}
