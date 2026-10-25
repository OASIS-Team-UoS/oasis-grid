// Extracted OASIS geometry: SPDX-License-Identifier: GPL-3.0-or-later
// Generated mesh-only wrapper; original function bodies are unchanged.
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
using real_t = double;
using realvector = std::vector<double>;
using intvector = std::vector<int>;
class Icogrid {
public:
 void sphere_ico(realvector&, int, int, int, int, int, int, intvector&, int);
 void find_qpoints_thor(const intvector&, const realvector&, const intvector&, realvector&) const;
 void relocate_centres_thor(const realvector&, const intvector&, realvector&) const;
};
namespace {
    inline real_t PI_const() {
        return static_cast<real_t>(3.141592653589793238462643383279502884L);
    }
    inline void normalize3(real_t& x, real_t& y, real_t& z) {
        const real_t n = std::sqrt(std::pow(x, static_cast<real_t>(2.0)) +
                                   std::pow(y, static_cast<real_t>(2.0)) +
                                   std::pow(z, static_cast<real_t>(2.0)));
        x /= n;
        y /= n;
        z /= n;
    }
}


void Icogrid::sphere_ico(realvector &xyz        , // Vertices coordinates (Cartesian; radius 1) [num*3]
                         int         glevel     , // Number of recursive iterations to increase horizontal resolution
                         int         n_region   , // Number of regions per main rhombus
                         int         nl_region  , // nl_region^2 is the number of points in the faces
                         int         nl2        , // nl2 is the number of points in the faces (redundant with nl_region^2)
                         int         kxl        , // kxl^2 is the number of faces inside the main rhombi
                         int         nfaces     , // Number of faces per base rhombus (should equal kxl^2)
                         intvector  &pent_ind   , // Pentagons' indexes (size >= 12)
                         int         num        ){// Total number of vertices (including the two poles)

    // --------------------------
    // Basic validation
    // --------------------------
    if (glevel < 0 || n_region <= 0 || nl_region <= 0 || nl2 <= 0 || kxl <= 0 || nfaces <= 0 || num <= 0) {
        throw std::runtime_error("sphere_ico: invalid non-positive inputs.");
    }
    // Check face/local consistency:
    //   nl2 must match nl_region^2
    //   n_region must match nfaces * nl2 (for divide_face > 0, n_region > nl2)
    {
        const int nl_region_sq = nl_region * nl_region;
        if (nl_region_sq != nl2 || n_region != nfaces * nl2) {
            std::ostringstream oss;
            oss << "sphere_ico: inconsistent n_region/nfaces/nl_region/nl2. "
                << "nl_region^2=" << nl_region_sq
                << " nfaces*nl2=" << (nfaces * nl2)
                << " n_region=" << n_region
                << " nl2=" << nl2;
            throw std::runtime_error(oss.str());
        }
    }
    // Check kxl^2 consistency
    {
        const int kxl_sq = kxl * kxl;
        if (kxl_sq != nfaces) {
            std::ostringstream oss;
            oss << "sphere_ico: inconsistent kxl/nfaces. kxl^2=" << kxl_sq << " nfaces=" << nfaces;
            throw std::runtime_error(oss.str());
        }
    }
    // Expect 10 base rhombi + 2 poles in pent_ind list (at least 12 entries)
    if (pent_ind.size() < 12u) {
        pent_ind.assign(12u, -1);
    }
    // xyz must have capacity for all points
    if (xyz.size() != static_cast<std::size_t>(num) * 3u) {
        xyz.assign(static_cast<std::size_t>(num) * 3u, static_cast<real_t>(0));
    }

    // --------------------------
    // Integer refinement sizes
    // --------------------------
    const int sizei   = static_cast<int>(std::pow(2.0, glevel)); // 2^glevel
    const int sizeip1 = sizei + 1;
    const int sizei2  = sizeip1 * sizeip1;

    // --------------------------
    // Temporary buffers
    // --------------------------
    intvector  rhombi(10 * sizei2, 0);              // rhombus (halo-inclusive) indexes into xyzi
    intvector  rhomb(10 * sizei * sizei, 0);        // rhombus (interior-only) indexes into xyzi
    realvector xyzi(static_cast<std::size_t>(3) * (10u * static_cast<std::size_t>(sizei2) + 2u),
                    static_cast<real_t>(0));        // temporary vertices (include two poles)

    // --------------------------
    // Icosahedron base geometry
    // --------------------------
    const real_t PI = PI_const();
    const real_t w  = static_cast<real_t>(2.0) * std::acos(static_cast<real_t>(1.0) /
                      (static_cast<real_t>(2.0) * std::sin(PI / static_cast<real_t>(5.0))));

    // Poles
    xyzi[0*3 + 0] = static_cast<real_t>(0);
    xyzi[0*3 + 1] = static_cast<real_t>(0);
    xyzi[0*3 + 2] = static_cast<real_t>(1);
    xyzi[1*3 + 0] = static_cast<real_t>(0);
    xyzi[1*3 + 1] = static_cast<real_t>(0);
    xyzi[1*3 + 2] = static_cast<real_t>(-1);

    // Other 10 points (indices 2..11)
    // 3
    xyzi[2*3 + 0] = std::cos(-PI/5)*std::cos(PI/2 - w);
    xyzi[2*3 + 1] = std::sin(-PI/5)*std::cos(PI/2 - w);
    xyzi[2*3 + 2] = std::sin(PI/2 - w);
    rhombi[0] = 2;
    // 4
    xyzi[3*3 + 0] = std::cos( PI/5)*std::cos(PI/2 - w);
    xyzi[3*3 + 1] = std::sin( PI/5)*std::cos(PI/2 - w);
    xyzi[3*3 + 2] = std::sin(PI/2 - w);
    rhombi[1] = 3;
    // 5
    xyzi[4*3 + 0] = std::cos( 3*PI/5)*std::cos(PI/2 - w);
    xyzi[4*3 + 1] = std::sin( 3*PI/5)*std::cos(PI/2 - w);
    xyzi[4*3 + 2] = std::sin(PI/2 - w);
    rhombi[2] = 4;
    // 6
    xyzi[5*3 + 0] = std::cos( PI)*std::cos(PI/2 - w);
    xyzi[5*3 + 1] = std::sin( PI)*std::cos(PI/2 - w);
    xyzi[5*3 + 2] = std::sin(PI/2 - w);
    rhombi[3] = 5;
    // 7
    xyzi[6*3 + 0] = std::cos(-(3.0/5.0)*PI)*std::cos(PI/2 - w);
    xyzi[6*3 + 1] = std::sin(-(3.0/5.0)*PI)*std::cos(PI/2 - w);
    xyzi[6*3 + 2] = std::sin(PI/2 - w);
    rhombi[4] = 6;
    // 8
    xyzi[7*3 + 0] = std::cos(0.0)*std::cos(w - PI/2);
    xyzi[7*3 + 1] = std::sin(0.0)*std::cos(w - PI/2);
    xyzi[7*3 + 2] = std::sin(w - PI/2);
    rhombi[5] = 7;
    // 9
    xyzi[8*3 + 0]  = std::cos(2*PI/5)*std::cos(w - PI/2);
    xyzi[8*3 + 1]  = std::sin(2*PI/5)*std::cos(w - PI/2);
    xyzi[8*3 + 2]  = std::sin(w - PI/2);
    rhombi[6] = 8;
    // 10
    xyzi[9*3 + 0]  = std::cos(4*PI/5)*std::cos(w - PI/2);
    xyzi[9*3 + 1]  = std::sin(4*PI/5)*std::cos(w - PI/2);
    xyzi[9*3 + 2]  = std::sin(w - PI/2);
    rhombi[7] = 9;
    // 11
    xyzi[10*3 + 0] = std::cos(-4*PI/5)*std::cos(w - PI/2);
    xyzi[10*3 + 1] = std::sin(-4*PI/5)*std::cos(w - PI/2);
    xyzi[10*3 + 2] = std::sin(w - PI/2);
    rhombi[8] = 10;
    // 12
    xyzi[11*3 + 0] = std::cos(-2*PI/5)*std::cos(w - PI/2);
    xyzi[11*3 + 1] = std::sin(-2*PI/5)*std::cos(w - PI/2);
    xyzi[11*3 + 2] = std::sin(w - PI/2);
    rhombi[9] = 11;

    // --------------------------
    // Standard grid point indexes for the 10 rhombi
    // --------------------------
    for (int faces = 0; faces < 10; faces++){
        if (faces == 0){
            rhombi[faces]                               = 2;
            rhombi[sizei*sizeip1*10 + faces]            = 0;
            rhombi[sizei*10 + faces]                    = 7;
            rhombi[sizei*sizeip1*10 + sizei*10 + faces] = 3;
        } else if (faces == 1){
            rhombi[faces]                               = 3;
            rhombi[sizei*sizeip1*10 + faces]            = 0;
            rhombi[sizei*10 + faces]                    = 8;
            rhombi[sizei*sizeip1*10 + sizei*10 + faces] = 4;
        } else if (faces == 2){
            rhombi[faces]                               = 4;
            rhombi[sizei*sizeip1*10 + faces]            = 0;
            rhombi[sizei*10 + faces]                    = 9;
            rhombi[sizei*sizeip1*10 + sizei*10 + faces] = 5;
        } else if (faces == 3){
            rhombi[faces]                               = 5;
            rhombi[sizei*sizeip1*10 + faces]            = 0;
            rhombi[sizei*10 + faces]                    = 10;
            rhombi[sizei*sizeip1*10 + sizei*10 + faces] = 6;
        } else if (faces == 4){
            rhombi[faces]                               = 6;
            rhombi[sizei*sizeip1*10 + faces]            = 0;
            rhombi[sizei*10 + faces]                    = 11;
            rhombi[sizei*sizeip1*10 + sizei*10 + faces] = 2;
        } else if (faces == 5){
            rhombi[faces]                               = 7;
            rhombi[sizei*sizeip1*10 + faces]            = 3;
            rhombi[sizei*10 + faces]                    = 1;
            rhombi[sizei*sizeip1*10 + sizei*10 + faces] = 8;
        } else if (faces == 6){
            rhombi[faces]                               = 8;
            rhombi[sizei*sizeip1*10 + faces]            = 4;
            rhombi[sizei*10 + faces]                    = 1;
            rhombi[sizei*sizeip1*10 + sizei*10 + faces] = 9;
        } else if (faces == 7){
            rhombi[faces]                               = 9;
            rhombi[sizei*sizeip1*10 + faces]            = 5;
            rhombi[sizei*10 + faces]                    = 1;
            rhombi[sizei*sizeip1*10 + sizei*10 + faces] = 10;
        } else if (faces == 8){
            rhombi[faces]                               = 10;
            rhombi[sizei*sizeip1*10 + faces]            = 6;
            rhombi[sizei*10 + faces]                    = 1;
            rhombi[sizei*sizeip1*10 + sizei*10 + faces] = 11;
        } else if (faces == 9){
            rhombi[faces]                               = 11;
            rhombi[sizei*sizeip1*10 + faces]            = 2;
            rhombi[sizei*10 + faces]                    = 1;
            rhombi[sizei*sizeip1*10 + sizei*10 + faces] = 7;
        }
    }

    // --------------------------
    // Recursive subdivision
    // --------------------------
    int count_points = 12; // initial icosahedron vertices

    // Edge midpoints (per rhombus)
    for (int faces = 0; faces < 10; faces++){
        for (int grd = 0; grd < glevel; grd++){
            const int ind_jump = static_cast<int>(std::pow(2.0, glevel - grd - 1));
            int ind = ind_jump;

            for (int it = 0; it < static_cast<int>(std::pow(2.0, grd)); it++){
                // Left edge (x varying)
                {
                    const int indx  = ind;
                    const int indy  = 0;
                    const int indx1 = ind + ind_jump;
                    const int indy1 = 0;
                    const int indx2 = ind - ind_jump;
                    const int indy2 = 0;
                    rhombi[indx*sizeip1*10 + indy*10 + faces] = count_points;
                    const int a = rhombi[indx1*sizeip1*10 + indy1*10 + faces];
                    const int b = rhombi[indx2*sizeip1*10 + indy2*10 + faces];
                    xyzi[count_points*3 + 0] = (xyzi[a*3 + 0] + xyzi[b*3 + 0]) * static_cast<real_t>(0.5);
                    xyzi[count_points*3 + 1] = (xyzi[a*3 + 1] + xyzi[b*3 + 1]) * static_cast<real_t>(0.5);
                    xyzi[count_points*3 + 2] = (xyzi[a*3 + 2] + xyzi[b*3 + 2]) * static_cast<real_t>(0.5);
                    normalize3(xyzi[count_points*3 + 0], xyzi[count_points*3 + 1], xyzi[count_points*3 + 2]);
                    ++count_points;
                }
                // Bottom edge (y varying)
                {
                    const int indx  = 0;
                    const int indy  = ind;
                    const int indx1 = 0;
                    const int indy1 = ind + ind_jump;
                    const int indx2 = 0;
                    const int indy2 = ind - ind_jump;
                    rhombi[indx*sizeip1*10 + indy*10 + faces] = count_points;
                    const int a = rhombi[indx1*sizeip1*10 + indy1*10 + faces];
                    const int b = rhombi[indx2*sizeip1*10 + indy2*10 + faces];
                    xyzi[count_points*3 + 0] = (xyzi[a*3 + 0] + xyzi[b*3 + 0]) * static_cast<real_t>(0.5);
                    xyzi[count_points*3 + 1] = (xyzi[a*3 + 1] + xyzi[b*3 + 1]) * static_cast<real_t>(0.5);
                    xyzi[count_points*3 + 2] = (xyzi[a*3 + 2] + xyzi[b*3 + 2]) * static_cast<real_t>(0.5);
                    normalize3(xyzi[count_points*3 + 0], xyzi[count_points*3 + 1], xyzi[count_points*3 + 2]);
                    ++count_points;
                }
                // Diagonal (center) of upper-left triangle
                {
                    const int indx  = ind;
                    const int indy  = ind;
                    const int indx1 = ind + ind_jump;
                    const int indy1 = ind + ind_jump;
                    const int indx2 = ind - ind_jump;
                    const int indy2 = ind - ind_jump;
                    rhombi[indx*sizeip1*10 + indy*10 + faces] = count_points;
                    const int a = rhombi[indx1*sizeip1*10 + indy1*10 + faces];
                    const int b = rhombi[indx2*sizeip1*10 + indy2*10 + faces];
                    xyzi[count_points*3 + 0] = (xyzi[a*3 + 0] + xyzi[b*3 + 0]) * static_cast<real_t>(0.5);
                    xyzi[count_points*3 + 1] = (xyzi[a*3 + 1] + xyzi[b*3 + 1]) * static_cast<real_t>(0.5);
                    xyzi[count_points*3 + 2] = (xyzi[a*3 + 2] + xyzi[b*3 + 2]) * static_cast<real_t>(0.5);
                    normalize3(xyzi[count_points*3 + 0], xyzi[count_points*3 + 1], xyzi[count_points*3 + 2]);
                    ++count_points;
                }
                // Right edge (x varying, y = sizei)
                {
                    const int indx  = ind;
                    const int indy  = sizei;
                    const int indx1 = ind + ind_jump;
                    const int indy1 = sizei;
                    const int indx2 = ind - ind_jump;
                    const int indy2 = sizei;
                    rhombi[indx*sizeip1*10 + indy*10 + faces] = count_points;
                    const int a = rhombi[indx1*sizeip1*10 + indy1*10 + faces];
                    const int b = rhombi[indx2*sizeip1*10 + indy2*10 + faces];
                    xyzi[count_points*3 + 0] = (xyzi[a*3 + 0] + xyzi[b*3 + 0]) * static_cast<real_t>(0.5);
                    xyzi[count_points*3 + 1] = (xyzi[a*3 + 1] + xyzi[b*3 + 1]) * static_cast<real_t>(0.5);
                    xyzi[count_points*3 + 2] = (xyzi[a*3 + 2] + xyzi[b*3 + 2]) * static_cast<real_t>(0.5);
                    normalize3(xyzi[count_points*3 + 0], xyzi[count_points*3 + 1], xyzi[count_points*3 + 2]);
                    ++count_points;
                }
                // Top edge (y varying, x = sizei)
                {
                    const int indx  = sizei;
                    const int indy  = ind;
                    const int indx1 = sizei;
                    const int indy1 = ind + ind_jump;
                    const int indx2 = sizei;
                    const int indy2 = ind - ind_jump;
                    rhombi[indx*sizeip1*10 + indy*10 + faces] = count_points;
                    const int a = rhombi[indx1*sizeip1*10 + indy1*10 + faces];
                    const int b = rhombi[indx2*sizeip1*10 + indy2*10 + faces];
                    xyzi[count_points*3 + 0] = (xyzi[a*3 + 0] + xyzi[b*3 + 0]) * static_cast<real_t>(0.5);
                    xyzi[count_points*3 + 1] = (xyzi[a*3 + 1] + xyzi[b*3 + 1]) * static_cast<real_t>(0.5);
                    xyzi[count_points*3 + 2] = (xyzi[a*3 + 2] + xyzi[b*3 + 2]) * static_cast<real_t>(0.5);
                    normalize3(xyzi[count_points*3 + 0], xyzi[count_points*3 + 1], xyzi[count_points*3 + 2]);
                    ++count_points;
                }

                ind += 2 * ind_jump;
            }
        }

        // Upper triangles (grd > 0)
        for (int grd = 0; grd < glevel; grd++){
            if (grd > 0){
                const int ind_jump = static_cast<int>(std::pow(2.0, glevel - grd - 1));
                int indj = ind_jump;
                for (int j = 0; j < (static_cast<int>(std::pow(2.0, grd)) - 1); j++){
                    int indy = ind_jump;
                    int indx = 0;
                    for (int i = 0; i < j + 1; i++){
                        const int rindx = sizei - indj + indx;
                        const int rindy = indy;

                        // center
                        rhombi[rindx*sizeip1*10 + rindy*10 + faces] = count_points;
                        {
                            const int a = rhombi[(rindx - ind_jump)*sizeip1*10 + (rindy - ind_jump)*10 + faces];
                            const int b = rhombi[(rindx + ind_jump)*sizeip1*10 + (rindy + ind_jump)*10 + faces];
                            xyzi[count_points*3 + 0] = (xyzi[a*3 + 0] + xyzi[b*3 + 0]) * static_cast<real_t>(0.5);
                            xyzi[count_points*3 + 1] = (xyzi[a*3 + 1] + xyzi[b*3 + 1]) * static_cast<real_t>(0.5);
                            xyzi[count_points*3 + 2] = (xyzi[a*3 + 2] + xyzi[b*3 + 2]) * static_cast<real_t>(0.5);
                            normalize3(xyzi[count_points*3 + 0], xyzi[count_points*3 + 1], xyzi[count_points*3 + 2]);
                            ++count_points;
                        }

                        // left
                        rhombi[(rindx - ind_jump)*sizeip1*10 + rindy*10 + faces] = count_points;
                        {
                            const int a = rhombi[(rindx - ind_jump)*sizeip1*10 + (rindy - ind_jump)*10 + faces];
                            const int b = rhombi[(rindx - ind_jump)*sizeip1*10 + (rindy + ind_jump)*10 + faces];
                            xyzi[count_points*3 + 0] = (xyzi[a*3 + 0] + xyzi[b*3 + 0]) * static_cast<real_t>(0.5);
                            xyzi[count_points*3 + 1] = (xyzi[a*3 + 1] + xyzi[b*3 + 1]) * static_cast<real_t>(0.5);
                            xyzi[count_points*3 + 2] = (xyzi[a*3 + 2] + xyzi[b*3 + 2]) * static_cast<real_t>(0.5);
                            normalize3(xyzi[count_points*3 + 0], xyzi[count_points*3 + 1], xyzi[count_points*3 + 2]);
                            ++count_points;
                        }

                        // top
                        rhombi[rindx*sizeip1*10 + (rindy + ind_jump)*10 + faces] = count_points;
                        {
                            const int a = rhombi[(rindx + ind_jump)*sizeip1*10 + (rindy + ind_jump)*10 + faces];
                            const int b = rhombi[(rindx - ind_jump)*sizeip1*10 + (rindy + ind_jump)*10 + faces];
                            xyzi[count_points*3 + 0] = (xyzi[a*3 + 0] + xyzi[b*3 + 0]) * static_cast<real_t>(0.5);
                            xyzi[count_points*3 + 1] = (xyzi[a*3 + 1] + xyzi[b*3 + 1]) * static_cast<real_t>(0.5);
                            xyzi[count_points*3 + 2] = (xyzi[a*3 + 2] + xyzi[b*3 + 2]) * static_cast<real_t>(0.5);
                            normalize3(xyzi[count_points*3 + 0], xyzi[count_points*3 + 1], xyzi[count_points*3 + 2]);
                            ++count_points;
                        }

                        indx += 2 * ind_jump;
                        indy += 2 * ind_jump;
                    }
                    indj += 2 * ind_jump;
                }
            }

            // Lower triangles (grd > 0)
            if (grd > 0){
                const int ind_jump = static_cast<int>(std::pow(2.0, glevel - grd - 1));
                int indj = ind_jump;
                for (int j = 0; j < (static_cast<int>(std::pow(2.0, grd)) - 1); j++){
                    int indy = 0;
                    int indx = ind_jump;
                    for (int i = 0; i < j + 1; i++){
                        const int rindx = indx;
                        const int rindy = sizei - indj + indy;

                        // center
                        rhombi[rindx*sizeip1*10 + rindy*10 + faces] = count_points;
                        {
                            const int a = rhombi[(rindx - ind_jump)*sizeip1*10 + (rindy - ind_jump)*10 + faces];
                            const int b = rhombi[(rindx + ind_jump)*sizeip1*10 + (rindy + ind_jump)*10 + faces];
                            xyzi[count_points*3 + 0] = (xyzi[a*3 + 0] + xyzi[b*3 + 0]) * static_cast<real_t>(0.5);
                            xyzi[count_points*3 + 1] = (xyzi[a*3 + 1] + xyzi[b*3 + 1]) * static_cast<real_t>(0.5);
                            xyzi[count_points*3 + 2] = (xyzi[a*3 + 2] + xyzi[b*3 + 2]) * static_cast<real_t>(0.5);
                            normalize3(xyzi[count_points*3 + 0], xyzi[count_points*3 + 1], xyzi[count_points*3 + 2]);
                            ++count_points;
                        }

                        // bottom
                        rhombi[rindx*sizeip1*10 + (rindy - ind_jump)*10 + faces] = count_points;
                        {
                            const int a = rhombi[(rindx - ind_jump)*sizeip1*10 + (rindy - ind_jump)*10 + faces];
                            const int b = rhombi[(rindx + ind_jump)*sizeip1*10 + (rindy - ind_jump)*10 + faces];
                            xyzi[count_points*3 + 0] = (xyzi[a*3 + 0] + xyzi[b*3 + 0]) * static_cast<real_t>(0.5);
                            xyzi[count_points*3 + 1] = (xyzi[a*3 + 1] + xyzi[b*3 + 1]) * static_cast<real_t>(0.5);
                            xyzi[count_points*3 + 2] = (xyzi[a*3 + 2] + xyzi[b*3 + 2]) * static_cast<real_t>(0.5);
                            normalize3(xyzi[count_points*3 + 0], xyzi[count_points*3 + 1], xyzi[count_points*3 + 2]);
                            ++count_points;
                        }

                        // right
                        rhombi[(rindx + ind_jump)*sizeip1*10 + rindy*10 + faces] = count_points;
                        {
                            const int a = rhombi[(rindx + ind_jump)*sizeip1*10 + (rindy + ind_jump)*10 + faces];
                            const int b = rhombi[(rindx + ind_jump)*sizeip1*10 + (rindy - ind_jump)*10 + faces];
                            xyzi[count_points*3 + 0] = (xyzi[a*3 + 0] + xyzi[b*3 + 0]) * static_cast<real_t>(0.5);
                            xyzi[count_points*3 + 1] = (xyzi[a*3 + 1] + xyzi[b*3 + 1]) * static_cast<real_t>(0.5);
                            xyzi[count_points*3 + 2] = (xyzi[a*3 + 2] + xyzi[b*3 + 2]) * static_cast<real_t>(0.5);
                            normalize3(xyzi[count_points*3 + 0], xyzi[count_points*3 + 1], xyzi[count_points*3 + 2]);
                            ++count_points;
                        }

                        indx += 2 * ind_jump;
                        indy += 2 * ind_jump;
                    }
                    indj += 2 * ind_jump;
                }
            }
        }
    }

    // --------------------------
    // Assemble interior-only rhomb indices
    // --------------------------
    const int nli_region = static_cast<int>(std::sqrt(n_region * 1.0));
    for (int fc = 0; fc < 10; fc++) {
        for (int j = 0; j < sizei; j++) {
            for (int i = 0; i < sizei; i++) {
                rhomb[i*sizei*10 + j*10 + fc] = rhombi[i*sizeip1*10 + j*10 + fc];
            }
        }
    }

    // --------------------------
    // Write final xyz (per-face tiled layout)
    // --------------------------
    for (int fc = 0; fc < 10; fc++){
        for (int kx = 0; kx < kxl; kx++){
            for (int ky = 0; ky < kxl; ky++){
                for (int i = 0; i < nl_region; i++){
                    for (int j = 0; j < nl_region; j++){
                        const int src_idx = rhomb[((ky*nl_region + j)*nli_region + kx*nl_region + i)*10 + fc];
                        const std::size_t dst = static_cast<std::size_t>(fc*nfaces*nl2 + ky*kxl*nl2 + kx*nl2 + j*nl_region + i) * 3u;
                        xyz[dst + 0] = xyzi[static_cast<std::size_t>(src_idx) * 3u + 0u];
                        xyz[dst + 1] = xyzi[static_cast<std::size_t>(src_idx) * 3u + 1u];
                        xyz[dst + 2] = xyzi[static_cast<std::size_t>(src_idx) * 3u + 2u];
                    }
                }
            }
        }
    }

    // --------------------------
    // Poles
    // --------------------------
    xyz[static_cast<std::size_t>(num - 2) * 3u + 0u] = static_cast<real_t>(0);
    xyz[static_cast<std::size_t>(num - 2) * 3u + 1u] = static_cast<real_t>(0);
    xyz[static_cast<std::size_t>(num - 2) * 3u + 2u] = static_cast<real_t>(1);

    xyz[static_cast<std::size_t>(num - 1) * 3u + 0u] = static_cast<real_t>(0);
    xyz[static_cast<std::size_t>(num - 1) * 3u + 1u] = static_cast<real_t>(0);
    xyz[static_cast<std::size_t>(num - 1) * 3u + 2u] = static_cast<real_t>(-1);

    // --------------------------
    // Pentagons' indexes (first cell of each face + poles)
    // --------------------------
    for (int faces = 0; faces < 10; faces++) {
        pent_ind[static_cast<std::size_t>(faces)] = faces * nfaces * nl2;
    }
    pent_ind[10] = num - 2; // north pole
    pent_ind[11] = num - 1; // south pole
}
void Icogrid::find_qpoints_thor(const intvector &point_local,
                                const realvector &xyz,
                                const intvector &pent_ind,
                                realvector &xyzq) const
{
    const int point_num = static_cast<int>(xyz.size() / 3u);
    if (point_num <= 0) return;
    if (point_local.size() < static_cast<std::size_t>(point_num) * 6u) return;

    xyzq.assign(static_cast<std::size_t>(point_num) * 6u * 3u, static_cast<real_t>(0));

    real_t v1[3];
    real_t v2[3];
    real_t vc1[3];
    real_t vc2[3];
    real_t vc3[3];
    real_t vgc[3];
    real_t l1, l2;
    int geo;

    for (int i = 0; i < point_num; i++){
        geo = 6; // Hexagons.
        for (int k = 0; k < 12; k++) if(i == pent_ind[static_cast<std::size_t>(k)]) geo = 5; // Pentagons.
        if(geo == 5){
            for (int j = 0; j < 4; j++){
                v2[0] = xyz[static_cast<std::size_t>(i)*3u + 0u];
                v2[1] = xyz[static_cast<std::size_t>(i)*3u + 1u];
                v2[2] = xyz[static_cast<std::size_t>(i)*3u + 2u];
                v1[0] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + static_cast<std::size_t>(j)])*3u + 0u];
                v1[1] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + static_cast<std::size_t>(j)])*3u + 1u];
                v1[2] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + static_cast<std::size_t>(j)])*3u + 2u];

                l1 = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];

                vc1[0] = v1[1]*v2[2] - v1[2]*v2[1];
                vc1[1] = v1[2]*v2[0] - v1[0]*v2[2];
                vc1[2] = v1[0]*v2[1] - v1[1]*v2[0];

                l2 = std::sqrt(vc1[0] * vc1[0] + vc1[1] * vc1[1] + vc1[2] * vc1[2]);

                vc1[0] = vc1[0] / l2 * std::atan2(l2, l1);
                vc1[1] = vc1[1] / l2 * std::atan2(l2, l1);
                vc1[2] = vc1[2] / l2 * std::atan2(l2, l1);

                v2[0] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + static_cast<std::size_t>(j)])*3u + 0u];
                v2[1] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + static_cast<std::size_t>(j)])*3u + 1u];
                v2[2] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + static_cast<std::size_t>(j)])*3u + 2u];
                v1[0] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + static_cast<std::size_t>(j + 1)])*3u + 0u];
                v1[1] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + static_cast<std::size_t>(j + 1)])*3u + 1u];
                v1[2] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + static_cast<std::size_t>(j + 1)])*3u + 2u];

                l1 = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];

                vc2[0] = v1[1]*v2[2] - v1[2]*v2[1];
                vc2[1] = v1[2]*v2[0] - v1[0]*v2[2];
                vc2[2] = v1[0]*v2[1] - v1[1]*v2[0];

                l2 = std::sqrt(vc2[0] * vc2[0] + vc2[1] * vc2[1] + vc2[2] * vc2[2]);

                vc2[0] = vc2[0] / l2 * std::atan2(l2, l1);
                vc2[1] = vc2[1] / l2 * std::atan2(l2, l1);
                vc2[2] = vc2[2] / l2 * std::atan2(l2, l1);

                v2[0] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + static_cast<std::size_t>(j + 1)])*3u + 0u];
                v2[1] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + static_cast<std::size_t>(j + 1)])*3u + 1u];
                v2[2] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + static_cast<std::size_t>(j + 1)])*3u + 2u];
                v1[0] = xyz[static_cast<std::size_t>(i)*3u + 0u];
                v1[1] = xyz[static_cast<std::size_t>(i)*3u + 1u];
                v1[2] = xyz[static_cast<std::size_t>(i)*3u + 2u];

                l1 = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];

                vc3[0] = v1[1]*v2[2] - v1[2]*v2[1];
                vc3[1] = v1[2]*v2[0] - v1[0]*v2[2];
                vc3[2] = v1[0]*v2[1] - v1[1]*v2[0];

                l2 = std::sqrt(vc3[0] * vc3[0] + vc3[1] * vc3[1] + vc3[2] * vc3[2]);

                vc3[0] = vc3[0] / l2 * std::atan2(l2, l1);
                vc3[1] = vc3[1] / l2 * std::atan2(l2, l1);
                vc3[2] = vc3[2] / l2 * std::atan2(l2, l1);

                vgc[0] = vc1[0] + vc2[0] + vc3[0];
                vgc[1] = vc1[1] + vc2[1] + vc3[1];
                vgc[2] = vc1[2] + vc2[2] + vc3[2];

                l2 = std::sqrt(vgc[0] * vgc[0] + vgc[1] * vgc[1] + vgc[2] * vgc[2]);

                const std::size_t out = (static_cast<std::size_t>(i) * 6u + static_cast<std::size_t>(j)) * 3u;
                xyzq[out + 0u] = vgc[0] / l2;
                xyzq[out + 1u] = vgc[1] / l2;
                xyzq[out + 2u] = vgc[2] / l2;
            }

            // j = 4 (neighbors 4 and 0)
            {
                v2[0] = xyz[static_cast<std::size_t>(i)*3u + 0u];
                v2[1] = xyz[static_cast<std::size_t>(i)*3u + 1u];
                v2[2] = xyz[static_cast<std::size_t>(i)*3u + 2u];
                v1[0] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + 4u])*3u + 0u];
                v1[1] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + 4u])*3u + 1u];
                v1[2] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + 4u])*3u + 2u];

                l1 = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];

                vc1[0] = v1[1]*v2[2] - v1[2]*v2[1];
                vc1[1] = v1[2]*v2[0] - v1[0]*v2[2];
                vc1[2] = v1[0]*v2[1] - v1[1]*v2[0];

                l2 = std::sqrt(vc1[0] * vc1[0] + vc1[1] * vc1[1] + vc1[2] * vc1[2]);

                vc1[0] = vc1[0] / l2 * std::atan2(l2, l1);
                vc1[1] = vc1[1] / l2 * std::atan2(l2, l1);
                vc1[2] = vc1[2] / l2 * std::atan2(l2, l1);

                v2[0] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + 4u])*3u + 0u];
                v2[1] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + 4u])*3u + 1u];
                v2[2] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + 4u])*3u + 2u];
                v1[0] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + 0u])*3u + 0u];
                v1[1] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + 0u])*3u + 1u];
                v1[2] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + 0u])*3u + 2u];

                l1 = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];

                vc2[0] = v1[1]*v2[2] - v1[2]*v2[1];
                vc2[1] = v1[2]*v2[0] - v1[0]*v2[2];
                vc2[2] = v1[0]*v2[1] - v1[1]*v2[0];

                l2 = std::sqrt(vc2[0] * vc2[0] + vc2[1] * vc2[1] + vc2[2] * vc2[2]);

                vc2[0] = vc2[0] / l2 * std::atan2(l2, l1);
                vc2[1] = vc2[1] / l2 * std::atan2(l2, l1);
                vc2[2] = vc2[2] / l2 * std::atan2(l2, l1);

                v2[0] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + 0u])*3u + 0u];
                v2[1] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + 0u])*3u + 1u];
                v2[2] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + 0u])*3u + 2u];
                v1[0] = xyz[static_cast<std::size_t>(i)*3u + 0u];
                v1[1] = xyz[static_cast<std::size_t>(i)*3u + 1u];
                v1[2] = xyz[static_cast<std::size_t>(i)*3u + 2u];

                l1 = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];

                vc3[0] = v1[1]*v2[2] - v1[2]*v2[1];
                vc3[1] = v1[2]*v2[0] - v1[0]*v2[2];
                vc3[2] = v1[0]*v2[1] - v1[1]*v2[0];

                l2 = std::sqrt(vc3[0] * vc3[0] + vc3[1] * vc3[1] + vc3[2] * vc3[2]);

                vc3[0] = vc3[0] / l2 * std::atan2(l2, l1);
                vc3[1] = vc3[1] / l2 * std::atan2(l2, l1);
                vc3[2] = vc3[2] / l2 * std::atan2(l2, l1);

                vgc[0] = vc1[0] + vc2[0] + vc3[0];
                vgc[1] = vc1[1] + vc2[1] + vc3[1];
                vgc[2] = vc1[2] + vc2[2] + vc3[2];

                l2 = std::sqrt(vgc[0] * vgc[0] + vgc[1] * vgc[1] + vgc[2] * vgc[2]);

                const std::size_t out = (static_cast<std::size_t>(i) * 6u + 4u) * 3u;
                xyzq[out + 0u] = vgc[0] / l2;
                xyzq[out + 1u] = vgc[1] / l2;
                xyzq[out + 2u] = vgc[2] / l2;
            }
        }
        else{
            for (int j = 0; j < 5; j++){
                v2[0] = xyz[static_cast<std::size_t>(i)*3u + 0u];
                v2[1] = xyz[static_cast<std::size_t>(i)*3u + 1u];
                v2[2] = xyz[static_cast<std::size_t>(i)*3u + 2u];
                v1[0] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + static_cast<std::size_t>(j)])*3u + 0u];
                v1[1] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + static_cast<std::size_t>(j)])*3u + 1u];
                v1[2] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + static_cast<std::size_t>(j)])*3u + 2u];

                l1 = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];

                vc1[0] = v1[1]*v2[2] - v1[2]*v2[1];
                vc1[1] = v1[2]*v2[0] - v1[0]*v2[2];
                vc1[2] = v1[0]*v2[1] - v1[1]*v2[0];

                l2 = std::sqrt(vc1[0] * vc1[0] + vc1[1] * vc1[1] + vc1[2] * vc1[2]);

                vc1[0] = vc1[0] / l2 * std::atan2(l2, l1);
                vc1[1] = vc1[1] / l2 * std::atan2(l2, l1);
                vc1[2] = vc1[2] / l2 * std::atan2(l2, l1);

                v2[0] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + static_cast<std::size_t>(j)])*3u + 0u];
                v2[1] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + static_cast<std::size_t>(j)])*3u + 1u];
                v2[2] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + static_cast<std::size_t>(j)])*3u + 2u];
                v1[0] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + static_cast<std::size_t>(j + 1)])*3u + 0u];
                v1[1] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + static_cast<std::size_t>(j + 1)])*3u + 1u];
                v1[2] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + static_cast<std::size_t>(j + 1)])*3u + 2u];

                l1 = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];

                vc2[0] = v1[1]*v2[2] - v1[2]*v2[1];
                vc2[1] = v1[2]*v2[0] - v1[0]*v2[2];
                vc2[2] = v1[0]*v2[1] - v1[1]*v2[0];

                l2 = std::sqrt(vc2[0] * vc2[0] + vc2[1] * vc2[1] + vc2[2] * vc2[2]);

                vc2[0] = vc2[0] / l2 * std::atan2(l2, l1);
                vc2[1] = vc2[1] / l2 * std::atan2(l2, l1);
                vc2[2] = vc2[2] / l2 * std::atan2(l2, l1);

                v2[0] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + static_cast<std::size_t>(j + 1)])*3u + 0u];
                v2[1] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + static_cast<std::size_t>(j + 1)])*3u + 1u];
                v2[2] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + static_cast<std::size_t>(j + 1)])*3u + 2u];
                v1[0] = xyz[static_cast<std::size_t>(i)*3u + 0u];
                v1[1] = xyz[static_cast<std::size_t>(i)*3u + 1u];
                v1[2] = xyz[static_cast<std::size_t>(i)*3u + 2u];

                l1 = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];

                vc3[0] = v1[1]*v2[2] - v1[2]*v2[1];
                vc3[1] = v1[2]*v2[0] - v1[0]*v2[2];
                vc3[2] = v1[0]*v2[1] - v1[1]*v2[0];

                l2 = std::sqrt(vc3[0] * vc3[0] + vc3[1] * vc3[1] + vc3[2] * vc3[2]);

                vc3[0] = vc3[0] / l2 * std::atan2(l2, l1);
                vc3[1] = vc3[1] / l2 * std::atan2(l2, l1);
                vc3[2] = vc3[2] / l2 * std::atan2(l2, l1);

                vgc[0] = vc1[0] + vc2[0] + vc3[0];
                vgc[1] = vc1[1] + vc2[1] + vc3[1];
                vgc[2] = vc1[2] + vc2[2] + vc3[2];

                l2 = std::sqrt(vgc[0] * vgc[0] + vgc[1] * vgc[1] + vgc[2] * vgc[2]);

                const std::size_t out = (static_cast<std::size_t>(i) * 6u + static_cast<std::size_t>(j)) * 3u;
                xyzq[out + 0u] = vgc[0] / l2;
                xyzq[out + 1u] = vgc[1] / l2;
                xyzq[out + 2u] = vgc[2] / l2;
            }

            // j = 5 (neighbors 5 and 0)
            {
                v2[0] = xyz[static_cast<std::size_t>(i)*3u + 0u];
                v2[1] = xyz[static_cast<std::size_t>(i)*3u + 1u];
                v2[2] = xyz[static_cast<std::size_t>(i)*3u + 2u];
                v1[0] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + 5u])*3u + 0u];
                v1[1] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + 5u])*3u + 1u];
                v1[2] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + 5u])*3u + 2u];

                l1 = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];

                vc1[0] = v1[1]*v2[2] - v1[2]*v2[1];
                vc1[1] = v1[2]*v2[0] - v1[0]*v2[2];
                vc1[2] = v1[0]*v2[1] - v1[1]*v2[0];

                l2 = std::sqrt(vc1[0] * vc1[0] + vc1[1] * vc1[1] + vc1[2] * vc1[2]);

                vc1[0] = vc1[0] / l2 * std::atan2(l2, l1);
                vc1[1] = vc1[1] / l2 * std::atan2(l2, l1);
                vc1[2] = vc1[2] / l2 * std::atan2(l2, l1);

                v2[0] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + 5u])*3u + 0u];
                v2[1] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + 5u])*3u + 1u];
                v2[2] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + 5u])*3u + 2u];
                v1[0] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + 0u])*3u + 0u];
                v1[1] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + 0u])*3u + 1u];
                v1[2] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + 0u])*3u + 2u];

                l1 = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];

                vc2[0] = v1[1]*v2[2] - v1[2]*v2[1];
                vc2[1] = v1[2]*v2[0] - v1[0]*v2[2];
                vc2[2] = v1[0]*v2[1] - v1[1]*v2[0];

                l2 = std::sqrt(vc2[0] * vc2[0] + vc2[1] * vc2[1] + vc2[2] * vc2[2]);

                vc2[0] = vc2[0] / l2 * std::atan2(l2, l1);
                vc2[1] = vc2[1] / l2 * std::atan2(l2, l1);
                vc2[2] = vc2[2] / l2 * std::atan2(l2, l1);

                v2[0] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + 0u])*3u + 0u];
                v2[1] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + 0u])*3u + 1u];
                v2[2] = xyz[static_cast<std::size_t>(point_local[static_cast<std::size_t>(i)*6u + 0u])*3u + 2u];
                v1[0] = xyz[static_cast<std::size_t>(i)*3u + 0u];
                v1[1] = xyz[static_cast<std::size_t>(i)*3u + 1u];
                v1[2] = xyz[static_cast<std::size_t>(i)*3u + 2u];

                l1 = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];

                vc3[0] = v1[1]*v2[2] - v1[2]*v2[1];
                vc3[1] = v1[2]*v2[0] - v1[0]*v2[2];
                vc3[2] = v1[0]*v2[1] - v1[1]*v2[0];

                l2 = std::sqrt(vc3[0] * vc3[0] + vc3[1] * vc3[1] + vc3[2] * vc3[2]);

                vc3[0] = vc3[0] / l2 * std::atan2(l2, l1);
                vc3[1] = vc3[1] / l2 * std::atan2(l2, l1);
                vc3[2] = vc3[2] / l2 * std::atan2(l2, l1);

                vgc[0] = vc1[0] + vc2[0] + vc3[0];
                vgc[1] = vc1[1] + vc2[1] + vc3[1];
                vgc[2] = vc1[2] + vc2[2] + vc3[2];

                l2 = std::sqrt(vgc[0] * vgc[0] + vgc[1] * vgc[1] + vgc[2] * vgc[2]);

                const std::size_t out = (static_cast<std::size_t>(i) * 6u + 5u) * 3u;
                xyzq[out + 0u] = vgc[0] / l2;
                xyzq[out + 1u] = vgc[1] / l2;
                xyzq[out + 2u] = vgc[2] / l2;
            }
        }
    }
}
void Icogrid::relocate_centres_thor(const realvector &xyzq,
                                    const intvector &pent_ind,
                                    realvector &xyz) const
{
    const int point_num = static_cast<int>(xyz.size() / 3u);
    if (point_num <= 0) return;
    if (xyzq.size() < static_cast<std::size_t>(point_num) * 6u * 3u) return;

    real_t v1[3];
    real_t v2[3];
    real_t vc1[3];
    real_t vc2[3];
    real_t vc3[3];
    real_t vc4[3];
    real_t vc5[3];
    real_t vc6[3];
    real_t vgc[3];
    real_t l1, l2;
    int geo;

    for (int i = 0; i < point_num; i++){
        geo = 6; // Hexagons.
        for (int k = 0; k < 12; k++) if(i == pent_ind[static_cast<std::size_t>(k)]) geo = 5; // Pentagons.
        if(geo == 5){

            v1[0] = xyzq[static_cast<std::size_t>(i)*6u*3u + 1u*3u + 0u];
            v1[1] = xyzq[static_cast<std::size_t>(i)*6u*3u + 1u*3u + 1u];
            v1[2] = xyzq[static_cast<std::size_t>(i)*6u*3u + 1u*3u + 2u];
            v2[0] = xyzq[static_cast<std::size_t>(i)*6u*3u + 0u*3u + 0u];
            v2[1] = xyzq[static_cast<std::size_t>(i)*6u*3u + 0u*3u + 1u];
            v2[2] = xyzq[static_cast<std::size_t>(i)*6u*3u + 0u*3u + 2u];

            l1 = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];

            vc1[0] = v1[1]*v2[2] - v1[2]*v2[1];
            vc1[1] = v1[2]*v2[0] - v1[0]*v2[2];
            vc1[2] = v1[0]*v2[1] - v1[1]*v2[0];

            l2 = std::sqrt(vc1[0] * vc1[0] + vc1[1] * vc1[1] + vc1[2] * vc1[2]);

            vc1[0] = vc1[0] / l2 * std::atan2(l2, l1);
            vc1[1] = vc1[1] / l2 * std::atan2(l2, l1);
            vc1[2] = vc1[2] / l2 * std::atan2(l2, l1);

            v1[0] = xyzq[static_cast<std::size_t>(i)*6u*3u + 2u*3u + 0u];
            v1[1] = xyzq[static_cast<std::size_t>(i)*6u*3u + 2u*3u + 1u];
            v1[2] = xyzq[static_cast<std::size_t>(i)*6u*3u + 2u*3u + 2u];
            v2[0] = xyzq[static_cast<std::size_t>(i)*6u*3u + 1u*3u + 0u];
            v2[1] = xyzq[static_cast<std::size_t>(i)*6u*3u + 1u*3u + 1u];
            v2[2] = xyzq[static_cast<std::size_t>(i)*6u*3u + 1u*3u + 2u];

            l1 = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];

            vc2[0] = v1[1]*v2[2] - v1[2]*v2[1];
            vc2[1] = v1[2]*v2[0] - v1[0]*v2[2];
            vc2[2] = v1[0]*v2[1] - v1[1]*v2[0];

            l2 = std::sqrt(vc2[0] * vc2[0] + vc2[1] * vc2[1] + vc2[2] * vc2[2]);

            vc2[0] = vc2[0] / l2 * std::atan2(l2, l1);
            vc2[1] = vc2[1] / l2 * std::atan2(l2, l1);
            vc2[2] = vc2[2] / l2 * std::atan2(l2, l1);

            v1[0] = xyzq[static_cast<std::size_t>(i)*6u*3u + 3u*3u + 0u];
            v1[1] = xyzq[static_cast<std::size_t>(i)*6u*3u + 3u*3u + 1u];
            v1[2] = xyzq[static_cast<std::size_t>(i)*6u*3u + 3u*3u + 2u];
            v2[0] = xyzq[static_cast<std::size_t>(i)*6u*3u + 2u*3u + 0u];
            v2[1] = xyzq[static_cast<std::size_t>(i)*6u*3u + 2u*3u + 1u];
            v2[2] = xyzq[static_cast<std::size_t>(i)*6u*3u + 2u*3u + 2u];

            l1 = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];

            vc3[0] = v1[1]*v2[2] - v1[2]*v2[1];
            vc3[1] = v1[2]*v2[0] - v1[0]*v2[2];
            vc3[2] = v1[0]*v2[1] - v1[1]*v2[0];

            l2 = std::sqrt(vc3[0] * vc3[0] + vc3[1] * vc3[1] + vc3[2] * vc3[2]);

            vc3[0] = vc3[0] / l2 * std::atan2(l2, l1);
            vc3[1] = vc3[1] / l2 * std::atan2(l2, l1);
            vc3[2] = vc3[2] / l2 * std::atan2(l2, l1);

            v1[0] = xyzq[static_cast<std::size_t>(i)*6u*3u + 4u*3u + 0u];
            v1[1] = xyzq[static_cast<std::size_t>(i)*6u*3u + 4u*3u + 1u];
            v1[2] = xyzq[static_cast<std::size_t>(i)*6u*3u + 4u*3u + 2u];
            v2[0] = xyzq[static_cast<std::size_t>(i)*6u*3u + 3u*3u + 0u];
            v2[1] = xyzq[static_cast<std::size_t>(i)*6u*3u + 3u*3u + 1u];
            v2[2] = xyzq[static_cast<std::size_t>(i)*6u*3u + 3u*3u + 2u];

            l1 = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];

            vc4[0] = v1[1]*v2[2] - v1[2]*v2[1];
            vc4[1] = v1[2]*v2[0] - v1[0]*v2[2];
            vc4[2] = v1[0]*v2[1] - v1[1]*v2[0];

            l2 = std::sqrt(vc4[0] * vc4[0] + vc4[1] * vc4[1] + vc4[2] * vc4[2]);

            vc4[0] = vc4[0] / l2 * std::atan2(l2, l1);
            vc4[1] = vc4[1] / l2 * std::atan2(l2, l1);
            vc4[2] = vc4[2] / l2 * std::atan2(l2, l1);

            v1[0] = xyzq[static_cast<std::size_t>(i)*6u*3u + 0u*3u + 0u];
            v1[1] = xyzq[static_cast<std::size_t>(i)*6u*3u + 0u*3u + 1u];
            v1[2] = xyzq[static_cast<std::size_t>(i)*6u*3u + 0u*3u + 2u];
            v2[0] = xyzq[static_cast<std::size_t>(i)*6u*3u + 4u*3u + 0u];
            v2[1] = xyzq[static_cast<std::size_t>(i)*6u*3u + 4u*3u + 1u];
            v2[2] = xyzq[static_cast<std::size_t>(i)*6u*3u + 4u*3u + 2u];

            l1 = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];

            vc5[0] = v1[1]*v2[2] - v1[2]*v2[1];
            vc5[1] = v1[2]*v2[0] - v1[0]*v2[2];
            vc5[2] = v1[0]*v2[1] - v1[1]*v2[0];

            l2 = std::sqrt(vc5[0] * vc5[0] + vc5[1] * vc5[1] + vc5[2] * vc5[2]);

            vc5[0] = vc5[0] / l2 * std::atan2(l2, l1);
            vc5[1] = vc5[1] / l2 * std::atan2(l2, l1);
            vc5[2] = vc5[2] / l2 * std::atan2(l2, l1);

            vgc[0] = vc1[0] + vc2[0] + vc3[0] + vc4[0] + vc5[0];
            vgc[1] = vc1[1] + vc2[1] + vc3[1] + vc4[1] + vc5[1];
            vgc[2] = vc1[2] + vc2[2] + vc3[2] + vc4[2] + vc5[2];

            l2 = std::sqrt(vgc[0] * vgc[0] + vgc[1] * vgc[1] + vgc[2] * vgc[2]);

            xyz[static_cast<std::size_t>(i)*3u + 0u] = vgc[0] / l2;
            xyz[static_cast<std::size_t>(i)*3u + 1u] = vgc[1] / l2;
            xyz[static_cast<std::size_t>(i)*3u + 2u] = vgc[2] / l2;
        }
        else{
            v1[0] = xyzq[static_cast<std::size_t>(i)*6u*3u + 1u*3u + 0u];
            v1[1] = xyzq[static_cast<std::size_t>(i)*6u*3u + 1u*3u + 1u];
            v1[2] = xyzq[static_cast<std::size_t>(i)*6u*3u + 1u*3u + 2u];
            v2[0] = xyzq[static_cast<std::size_t>(i)*6u*3u + 0u*3u + 0u];
            v2[1] = xyzq[static_cast<std::size_t>(i)*6u*3u + 0u*3u + 1u];
            v2[2] = xyzq[static_cast<std::size_t>(i)*6u*3u + 0u*3u + 2u];

            l1 = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];

            vc1[0] = v1[1]*v2[2] - v1[2]*v2[1];
            vc1[1] = v1[2]*v2[0] - v1[0]*v2[2];
            vc1[2] = v1[0]*v2[1] - v1[1]*v2[0];

            l2 = std::sqrt(vc1[0] * vc1[0] + vc1[1] * vc1[1] + vc1[2] * vc1[2]);

            vc1[0] = vc1[0] / l2 * std::atan2(l2, l1);
            vc1[1] = vc1[1] / l2 * std::atan2(l2, l1);
            vc1[2] = vc1[2] / l2 * std::atan2(l2, l1);

            v1[0] = xyzq[static_cast<std::size_t>(i)*6u*3u + 2u*3u + 0u];
            v1[1] = xyzq[static_cast<std::size_t>(i)*6u*3u + 2u*3u + 1u];
            v1[2] = xyzq[static_cast<std::size_t>(i)*6u*3u + 2u*3u + 2u];
            v2[0] = xyzq[static_cast<std::size_t>(i)*6u*3u + 1u*3u + 0u];
            v2[1] = xyzq[static_cast<std::size_t>(i)*6u*3u + 1u*3u + 1u];
            v2[2] = xyzq[static_cast<std::size_t>(i)*6u*3u + 1u*3u + 2u];

            l1 = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];

            vc2[0] = v1[1]*v2[2] - v1[2]*v2[1];
            vc2[1] = v1[2]*v2[0] - v1[0]*v2[2];
            vc2[2] = v1[0]*v2[1] - v1[1]*v2[0];

            l2 = std::sqrt(vc2[0] * vc2[0] + vc2[1] * vc2[1] + vc2[2] * vc2[2]);

            vc2[0] = vc2[0] / l2 * std::atan2(l2, l1);
            vc2[1] = vc2[1] / l2 * std::atan2(l2, l1);
            vc2[2] = vc2[2] / l2 * std::atan2(l2, l1);

            v1[0] = xyzq[static_cast<std::size_t>(i)*6u*3u + 3u*3u + 0u];
            v1[1] = xyzq[static_cast<std::size_t>(i)*6u*3u + 3u*3u + 1u];
            v1[2] = xyzq[static_cast<std::size_t>(i)*6u*3u + 3u*3u + 2u];
            v2[0] = xyzq[static_cast<std::size_t>(i)*6u*3u + 2u*3u + 0u];
            v2[1] = xyzq[static_cast<std::size_t>(i)*6u*3u + 2u*3u + 1u];
            v2[2] = xyzq[static_cast<std::size_t>(i)*6u*3u + 2u*3u + 2u];

            l1 = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];

            vc3[0] = v1[1]*v2[2] - v1[2]*v2[1];
            vc3[1] = v1[2]*v2[0] - v1[0]*v2[2];
            vc3[2] = v1[0]*v2[1] - v1[1]*v2[0];

            l2 = std::sqrt(vc3[0] * vc3[0] + vc3[1] * vc3[1] + vc3[2] * vc3[2]);

            vc3[0] = vc3[0] / l2 * std::atan2(l2, l1);
            vc3[1] = vc3[1] / l2 * std::atan2(l2, l1);
            vc3[2] = vc3[2] / l2 * std::atan2(l2, l1);

            v1[0] = xyzq[static_cast<std::size_t>(i)*6u*3u + 4u*3u + 0u];
            v1[1] = xyzq[static_cast<std::size_t>(i)*6u*3u + 4u*3u + 1u];
            v1[2] = xyzq[static_cast<std::size_t>(i)*6u*3u + 4u*3u + 2u];
            v2[0] = xyzq[static_cast<std::size_t>(i)*6u*3u + 3u*3u + 0u];
            v2[1] = xyzq[static_cast<std::size_t>(i)*6u*3u + 3u*3u + 1u];
            v2[2] = xyzq[static_cast<std::size_t>(i)*6u*3u + 3u*3u + 2u];

            l1 = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];

            vc4[0] = v1[1]*v2[2] - v1[2]*v2[1];
            vc4[1] = v1[2]*v2[0] - v1[0]*v2[2];
            vc4[2] = v1[0]*v2[1] - v1[1]*v2[0];

            l2 = std::sqrt(vc4[0] * vc4[0] + vc4[1] * vc4[1] + vc4[2] * vc4[2]);

            vc4[0] = vc4[0] / l2 * std::atan2(l2, l1);
            vc4[1] = vc4[1] / l2 * std::atan2(l2, l1);
            vc4[2] = vc4[2] / l2 * std::atan2(l2, l1);

            v1[0] = xyzq[static_cast<std::size_t>(i)*6u*3u + 5u*3u + 0u];
            v1[1] = xyzq[static_cast<std::size_t>(i)*6u*3u + 5u*3u + 1u];
            v1[2] = xyzq[static_cast<std::size_t>(i)*6u*3u + 5u*3u + 2u];
            v2[0] = xyzq[static_cast<std::size_t>(i)*6u*3u + 4u*3u + 0u];
            v2[1] = xyzq[static_cast<std::size_t>(i)*6u*3u + 4u*3u + 1u];
            v2[2] = xyzq[static_cast<std::size_t>(i)*6u*3u + 4u*3u + 2u];

            l1 = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];

            vc5[0] = v1[1]*v2[2] - v1[2]*v2[1];
            vc5[1] = v1[2]*v2[0] - v1[0]*v2[2];
            vc5[2] = v1[0]*v2[1] - v1[1]*v2[0];

            l2 = std::sqrt(vc5[0] * vc5[0] + vc5[1] * vc5[1] + vc5[2] * vc5[2]);

            vc5[0] = vc5[0] / l2 * std::atan2(l2, l1);
            vc5[1] = vc5[1] / l2 * std::atan2(l2, l1);
            vc5[2] = vc5[2] / l2 * std::atan2(l2, l1);

            v1[0] = xyzq[static_cast<std::size_t>(i)*6u*3u + 0u*3u + 0u];
            v1[1] = xyzq[static_cast<std::size_t>(i)*6u*3u + 0u*3u + 1u];
            v1[2] = xyzq[static_cast<std::size_t>(i)*6u*3u + 0u*3u + 2u];
            v2[0] = xyzq[static_cast<std::size_t>(i)*6u*3u + 5u*3u + 0u];
            v2[1] = xyzq[static_cast<std::size_t>(i)*6u*3u + 5u*3u + 1u];
            v2[2] = xyzq[static_cast<std::size_t>(i)*6u*3u + 5u*3u + 2u];

            l1 = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];

            vc6[0] = v1[1]*v2[2] - v1[2]*v2[1];
            vc6[1] = v1[2]*v2[0] - v1[0]*v2[2];
            vc6[2] = v1[0]*v2[1] - v1[1]*v2[0];

            l2 = std::sqrt(vc6[0] * vc6[0] + vc6[1] * vc6[1] + vc6[2] * vc6[2]);

            vc6[0] = vc6[0] / l2 * std::atan2(l2, l1);
            vc6[1] = vc6[1] / l2 * std::atan2(l2, l1);
            vc6[2] = vc6[2] / l2 * std::atan2(l2, l1);

            vgc[0] = vc1[0] + vc2[0] + vc3[0] + vc4[0] + vc5[0] + vc6[0];
            vgc[1] = vc1[1] + vc2[1] + vc3[1] + vc4[1] + vc5[1] + vc6[1];
            vgc[2] = vc1[2] + vc2[2] + vc3[2] + vc4[2] + vc5[2] + vc6[2];

            l2 = std::sqrt(vgc[0] * vgc[0] + vgc[1] * vgc[1] + vgc[2] * vgc[2]);

            xyz[static_cast<std::size_t>(i)*3u + 0u] = vgc[0] / l2;
            xyz[static_cast<std::size_t>(i)*3u + 1u] = vgc[1] / l2;
            xyz[static_cast<std::size_t>(i)*3u + 2u] = vgc[2] / l2;
        }
    }
}

int main(int argc, char** argv) {
 if (argc != 3) return 2;
 const int level = std::stoi(argv[1]);
 const int n = 1 << level, count = 10*n*n+2;
 Icogrid grid;
 realvector xyz(3*count), q;
 intvector pent(12), rings(6*count, -1);
 grid.sphere_ico(xyz, level, n*n, n, n*n, 1, 1, pent, count);
 std::cout << std::setprecision(17) << std::scientific;
 if (std::string(argv[2]) == "centres") {
   for (int i=0; i<count; ++i)
     std::cout << xyz[3*i] << ' ' << xyz[3*i+1] << ' ' << xyz[3*i+2] << '\n';
   return 0;
 }
 for (auto& v : rings) if (!(std::cin >> v)) return 3;
 grid.find_qpoints_thor(rings, xyz, pent, q);
 realvector seed(xyz);
 grid.relocate_centres_thor(q, pent, xyz);
 for (int i=0; i<count; ++i) {
   for (int k=0; k<3; ++k) std::cout << seed[3*i+k] << ' ';
   for (int k=0; k<3; ++k) std::cout << xyz[3*i+k] << ' ';
   for (int k=0; k<18; ++k) std::cout << q[18*i+k] << ' ';
   std::cout << '\n';
 }
}
