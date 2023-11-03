#include "TMath.h"

TLorentzVector getNeutrinoP4(double MetPt, double MetPhi, TLorentzVector lep, int lepType){
        float MW_=80.385;

        double leppt = lep.Pt();
        double lepphi = lep.Phi();
        double lepeta = lep.Eta();
        double lepenergy = lep.Energy();

        double metpt = MetPt;
        double metphi = MetPhi;

        double  px = metpt*cos(metphi);
        double  py = metpt*sin(metphi);
        double  pz = 0;
        double  pxl= leppt*cos(lepphi);
        double  pyl= leppt*sin(lepphi);
        double  pzl= leppt*sinh(lepeta);
        double  El = lepenergy;
        double  a = pow(MW_,2) + pow(px+pxl,2) + pow(py+pyl,2) - px*px - py*py - El*El + pzl*pzl;
        double  b = 2.*pzl;
        double  A = b*b -4.*El*El;
        double  B = 2.*a*b;
        double  C = a*a-4.*(px*px+py*py)*El*El;

        ///////////////////////////pz for fnal
        double M_mu =  0;

        //if(lepType==1)M_mu=0.105658367;//mu
        //if(lepType==0)M_mu=0.00051099891;//electron

        int type=2; // use the small abs real root

        a = MW_*MW_ - M_mu*M_mu + 2.0*pxl*px + 2.0*pyl*py;
        A = 4.0*(El*El - pzl*pzl);
        B = -4.0*a*pzl;
        C = 4.0*El*El*(px*px + py*py) - a*a;

        double tmproot = B*B - 4.0*A*C;

        if (tmproot<0) {
            //std::cout << "Complex root detected, taking real part..." << std::endl;
            pz = - B/(2*A); // take real part of complex roots
        }
        else {
            double tmpsol1 = (-B + sqrt(tmproot))/(2.0*A);
            double tmpsol2 = (-B - sqrt(tmproot))/(2.0*A);
            //std::cout << " Neutrino Solutions: " << tmpsol1 << ", " << tmpsol2 << std::endl;

            if (type == 0 ) {
                // two real roots, pick the one closest to pz of muon
                if (TMath::Abs(tmpsol2-pzl) < TMath::Abs(tmpsol1-pzl)) { pz = tmpsol2; }
                else { pz = tmpsol1; }
                // if pz is > 300 pick the most central root
                if ( abs(pz) > 300. ) {
                    if (TMath::Abs(tmpsol1)<TMath::Abs(tmpsol2) ) { pz = tmpsol1; }
                    else { pz = tmpsol2; }
                }
            }
            if (type == 1 ) {
                // two real roots, pick the one closest to pz of muon
                if (TMath::Abs(tmpsol2-pzl) < TMath::Abs(tmpsol1-pzl)) { pz = tmpsol2; }
                else {pz = tmpsol1; }
            }
            if (type == 2 ) {
                // pick the most central root.
                if (TMath::Abs(tmpsol1)<TMath::Abs(tmpsol2) ) { pz = tmpsol1; }
                else { pz = tmpsol2; }
            }
            /*if (type == 3 ) {
             // pick the largest value of the cosine
             TVector3 p3w, p3mu;
             p3w.SetXYZ(pxl+px, pyl+py, pzl+ tmpsol1);
             p3mu.SetXYZ(pxl, pyl, pzl );

             double sinthcm1 = 2.*(p3mu.Perp(p3w))/MW_;
             p3w.SetXYZ(pxl+px, pyl+py, pzl+ tmpsol2);
             double sinthcm2 = 2.*(p3mu.Perp(p3w))/MW_;

             double costhcm1 = sqrt(1. - sinthcm1*sinthcm1);
             double costhcm2 = sqrt(1. - sinthcm2*sinthcm2);

             if ( costhcm1 > costhcm2 ) { pz = tmpsol1; otherSol_ = tmpsol2; }
             else { pz = tmpsol2;otherSol_ = tmpsol1; }

             }*///end of type3

        }//endl of if real root

        //dont correct pt neutrino
        TLorentzVector outP4;
        outP4.SetPxPyPzE(px,py,pz,sqrt(px*px+py*py+pz*pz));
        return outP4;

}