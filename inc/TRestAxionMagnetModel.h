/*************************************************************************
 * This file is part of the REST software framework.                     *
 *                                                                       *
 * Copyright (C) 2016 GIFNA/TREX (University of Zaragoza)                *
 * For more information see https://gifna.unizar.es/trex                 *
 *                                                                       *
 * REST is free software: you can redistribute it and/or modify          *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation, either version 3 of the License, or     *
 * (at your option) any later version.                                   *
 *                                                                       *
 * REST is distributed in the hope that it will be useful,               *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          *
 * GNU General Public License for more details.                          *
 *                                                                       *
 * You should have a copy of the GNU General Public License along with   *
 * REST in $REST_PATH/LICENSE.                                           *
 * If not, see https://www.gnu.org/licenses/.                            *
 * For the list of contributors see $REST_PATH/CREDITS.                  *
 *************************************************************************/

#ifndef REST_TRestAxionMagnetModel
#define REST_TRestAxionMagnetModel

#include <TRestAxionMagneticFit.h>
#include <TRestMetadata.h>

/// A class to fit and store the fit parameters from the magnetic field map
class TRestAxionMagnetModel : public TObject {
   private:
    std::vector<TRestAxionMagneticFit> fLinearFit;

   public:
    virtual std::vector<Double_t> GetComponentAlongPath(Int_t axis, TVector3 from, TVector3 to,
                                                        Double_t dl = 1., Int_t Nmax = 0) = 0;

    void Test(Double_t X = 340, Double_t Y = 0, Double_t dl = 50);

    const TRestAxionMagneticFit& GetFit(size_t n) const { return fLinearFit[n]; }

    ClassDef(TRestAxionMagnetModel, 1);
};
#endif
