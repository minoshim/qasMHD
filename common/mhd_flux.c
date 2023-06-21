#include "mhd_flux.h"
#include <math.h>
#include "common_func.h"

#define FL_EPS (1e-8)

void calc_flux_mlau(double rol, double vnl, double vtl, double vul, double btl, double bul, double prl,
		    double ror, double vnr, double vtr, double vur, double btr, double bur, double prr,
		    double bnc, double gamma, const double *dv,
		    double *fro, double *fmn, double *fmt, double *fmu, double *fbt, double *fbu, double *fen)
/* Calculate MLAU fluxes */
/* rol,vnl,vtl,vul,btl,bul,prl: input primitive variables at the left side */
/* ror,vnr,vtr,vur,btr,bur,prr: input primitive variables at the right side */
/* bnc: input normal magnetic field at the interface */
/* gamma: specific heat ratio */
/* fro,fmn,fmt,fmu,fbt,fbu,fen: output MLAU fluxes at the interface*/

/* Shock detection using dv */
{
  double gammam1i=1.0/(gamma-1.0);
  /* Bn at the interface */
  double bnc2=bnc*bnc;
  double sgn=(bnc > 0)?(1.0):(-1.0);
  double bnca=sgn*bnc;
  /* Variables at the left-face */
  double roli=1.0/rol;
  double vl2=vnl*vnl+vtl*vtl+vul*vul;
  double pml=0.5*(btl*btl+bul*bul);
  double ptl=prl+pml;
  double enl=gammam1i*prl+pml+0.5*rol*vl2;
  double hnl=(enl+ptl)*roli;
  /* Variables at the right-face */
  double rori=1.0/ror;
  double vr2=vnr*vnr+vtr*vtr+vur*vur;
  double pmr=0.5*(btr*btr+bur*bur);
  double ptr=prr+pmr;
  double enr=gammam1i*prr+pmr+0.5*ror*vr2;
  double hnr=(enr+ptr)*rori;
  /* Wave speeds */
  double cl2=gamma*prl*roli;
  double cr2=gamma*prr*rori;
  double cal2=bnc2*roli;
  double car2=bnc2*rori;
  double cbl2=cl2+cal2+2.0*pml*roli;
  double cbr2=cr2+car2+2.0*pmr*rori;
  double cfl2=0.5*(cbl2+sqrt(fabs(cbl2*cbl2-4.0*cl2*cal2)));
  double cfr2=0.5*(cbr2+sqrt(fabs(cbr2*cbr2-4.0*cr2*car2)));
  cbl2=vl2+cal2+2.0*pml*roli; /* Sound vel => convective vel. */
  cbr2=vr2+car2+2.0*pmr*rori; /* Sound vel => convective vel. */
  double ccl2=0.5*(cbl2+sqrt(fabs(cbl2*cbl2-4.0*vl2*cal2))); /* Sound vel => convective vel. */
  double ccr2=0.5*(cbr2+sqrt(fabs(cbr2*cbr2-4.0*vr2*car2))); /* Sound vel => convective vel. */
  double cmax=sqrt(max(cfl2,cfr2));
  double cfn=cmax;
  double cfni=1.0/cfn;
  double sl=min(0.0,min(vnl,vnr)-cmax);
  double sr=max(0.0,max(vnl,vnr)+cmax);

  /* Fast mach number functions */
  double mal=vnl*cfni;
  double mar=vnr*cfni;
  double snl=(mal > 0)?1.0:(-1.0);
  double snr=(mar > 0)?1.0:(-1.0);
  double theta=min(1.0,(cfn-min(dv[0],0.0))/(cfn-min(dv[1],0.0))); /* Shock detection */
  theta*=theta;
  theta*=theta;
  double alpha=0.1875;
  double beta=0.125;
  double malp1=mal+1.0;
  double marm1=mar-1.0;
  double mal2m1=mal*mal-1.0;
  double mar2m1=mar*mar-1.0;
  double bfl,ffl,bfr,ffr;
  if (snl*mal > 1.0){
    bfl=0.5*(1.0+snl);
    ffl=bfl*mal;
  } else{
    bfl=+0.25*malp1*malp1*(2.0-mal)+alpha*mal*mal2m1*mal2m1;
    ffl=+0.25*malp1*malp1+beta*mal2m1*mal2m1;
  }
  if (snr*mar > 1.0){
    bfr=0.5*(1.0-snr);
    ffr=bfr*mar;
  } else{
    bfr=+0.25*marm1*marm1*(2.0+mar)-alpha*mar*mar2m1*mar2m1;
    ffr=-0.25*marm1*marm1-beta*mar2m1*mar2m1;
  }

  /* Mass and Pressure fluxes */
  double mnc=ffl+ffr;
  double mflux=mnc*cfn;
  mflux-=max(1.0-fabs(mnc),0.0)*theta*cfni*(ptr-ptl)/(rol+ror);
  double dir=(mflux > 0)?1.0:(-1.0);
  double dirl=0.5*(1.0+dir);
  double dirr=1.0-dirl;
  mflux*=(dirl*rol+dirr*ror);
  double mfl=mflux*dirl;
  double mfr=mflux-mfl;
  double ptot=bfl*ptl+bfr*ptr;	/* AUSM(raw) */
  double cfmod=sqrt(max(ccl2,ccr2));		 /* Factor for all-speed */
  ptot-=0.25*bfl*bfr*(rol+ror)*cfmod*(vnr-vnl); /* AUSM+-up-like correction */
  /* ptot-=0.5*(1.0-min(1.0,cfmod*cfni))*(bfl+bfr-1.0)*(ptl+ptr); /\* SLAU-like correction *\/ */
  ptot-=0.5*(1.0-cfmod*cfni)*(bfl+bfr-1.0)*(ptl+ptr); /* SLAU2-like correction */

  /* Normal velocity */
  double slvl=sl-vnl;
  double srvr=sr-vnr;
  double snc=dirl*sl+dirr*sr;
  double vndnm=mflux+dirl*rol*slvl+dirr*ror*srvr;
  double vnc=((fabs(vndnm) > FL_EPS) && (snc != 0))?(mflux*snc/vndnm):(dirl*vnl+dirr*vnr);

  /* HLLD solution for tangential components in the Riemann fan */
  double slvc=sl-vnc;
  double srvc=sr-vnc;
  double slvci=1.0/slvc;
  double srvci=1.0/srvc;
  double rhdl=rol*slvl*slvc-bnc2;
  double rhdr=ror*srvr*srvc-bnc2;
  double dvtl=0.0,dvul=0.0,dbtl=0.0,dbul=0.0;
  double dvtr=0.0,dvur=0.0,dbtr=0.0,dbur=0.0;
  if (fabs(rhdl) > FL_EPS){
    double dvrhdli=(vnc-vnl)/rhdl;
    dvtl=-bnc*btl*dvrhdli;
    dvul=-bnc*bul*dvrhdli;
    dbtl=bnc2*btl*dvrhdli*slvci;
    dbul=bnc2*bul*dvrhdli*slvci;
  }
  if (fabs(rhdr) > FL_EPS){
    double dvrhdri=(vnc-vnr)/rhdr;
    dvtr=-bnc*btr*dvrhdri;
    dvur=-bnc*bur*dvrhdri;
    dbtr=bnc2*btr*dvrhdri*srvci;
    dbur=bnc2*bur*dvrhdri*srvci;
  }
  double cmpl=slvl*slvci;
  double cmpr=srvr*srvci;
  double ro2l=rol*cmpl;
  double srol=sqrt(ro2l);
  double vt2l=dvtl+vtl;
  double vu2l=dvul+vul;
  double bt2l=dbtl+btl*cmpl;
  double bu2l=dbul+bul*cmpl;
  double ro2r=ror*cmpr;
  double sror=sqrt(ro2r);
  double vt2r=dvtr+vtr;
  double vu2r=dvur+vur;
  double bt2r=dbtr+btr*cmpr;
  double bu2r=dbur+bur*cmpr;

  /* Magnetic tension flux */
  double sroi=1.0/(srol+sror);
  double srollri=srol*sroi;
  double srorlri=sror*sroi;
  double bnflb=sgn*min(max(0.0,srollri*(bnca+sror*vnc)),bnca);
  double bnfrb=sgn*min(max(0.0,srorlri*(bnca-srol*vnc)),bnca);
  double bnflv=sgn*min(max(0.0,srorlri*(bnca+mflux/sror)),bnca);
  double bnfrv=sgn*min(max(0.0,srollri*(bnca-mflux/srol)),bnca);
  double bnsrv=max(bnca-dir*vnc*(dirl*srol+dirr*sror),0.0);
  double ccb=sroi*bnsrv;
  double ccv=srol*sror*ccb;
  double bnvt=bnflb*vt2l+bnfrb*vt2r+ccb*(bt2r-bt2l);
  double bnvu=bnflb*vu2l+bnfrb*vu2r+ccb*(bu2r-bu2l);
  double bnbt=bnflv*bt2l+bnfrv*bt2r+ccv*(vt2r-vt2l);
  double bnbu=bnflv*bu2l+bnfrv*bu2r+ccv*(vu2r-vu2l);
  bnvt-=  vnc*(dirl*dbtl+dirr*dbtr);
  bnvu-=  vnc*(dirl*dbul+dirr*dbur);
  bnbt-=mflux*(dirl*dvtl+dirr*dvtr);
  bnbu-=mflux*(dirl*dvul+dirr*dvur);
  /* double bnvb=(bnc == 0)?0:(bnbt*bnvt+bnbu*bnvu)/bnc; */

  /* AUSM fluxes */
  *fro=mflux;
  *fmn=mfl*vnl     +mfr*vnr     +ptot-0.5*bnc2;
  *fmt=mfl*vtl     +mfr*vtr     -bnbt;
  *fmu=mfl*vul     +mfr*vur     -bnbu;
  *fbt=mfl*btl*roli+mfr*btr*rori-bnvt;
  *fbu=mfl*bul*roli+mfr*bur*rori-bnvu;

  /* Poyning flux term consistent with HLLD */
  double vb1=dirl*( vtl*btl + vul*bul) +dirr*( vtr*btr + vur*bur) ;
  double vb2=dirl*(vt2l*bt2l+vu2l*bu2l)+dirr*(vt2r*bt2r+vu2r*bu2r);
  double rv2b2=(mflux*vnc-bnc2);
  rv2b2*=rv2b2;
  double vb3=0.0;
  if (rv2b2 > FL_EPS*FL_EPS){
    vb3=(+((*fmt)*vnc+(*fbt)*bnc)*((*fmt)*bnc+(*fbt)*mflux)
	 +((*fmu)*vnc+(*fbu)*bnc)*((*fmu)*bnc+(*fbu)*mflux))/rv2b2;
  }
  double bnvb=sgn*(bnca*(snc*vb2-vnc*vb1)*(dirl*slvci+dirr*srvci)+bnsrv*(vb3-vb2));
  
  /* AUSM fluxes for energy */
  *fen=mfl*hnl     +mfr*hnr     -bnvb;
}

void calc_flux_roe(double rol, double vnl, double vtl, double vul, double btl, double bul, double prl,
		   double ror, double vnr, double vtr, double vur, double btr, double bur, double prr,
		   double bnc, double gamma, const double *dv,
		   double *fro, double *fmn, double *fmt, double *fmu, double *fbt, double *fbu, double *fen)
/* Calculate Roe fluxes */
/* rol,vnl,vtl,vul,btl,bul,prl: input primitive variables at the left side */
/* ror,vnr,vtr,vur,btr,bur,prr: input primitive variables at the right side */
/* bnc: input normal magnetic field at the interface */
/* gamma: specific heat ratio */
/* fro,fmn,fmt,fmu,fbt,fbu,fen: output Roe fluxes at the interface*/

/* dv is dummy (not used) */
{
  /* Adiabatic index */
  double gammam1=gamma-1.0;
  double gammam2=gamma-2.0;
  double gammam1i=1.0/gammam1;
  /* Bn at the interface */
  double bnc2=bnc*bnc;
  /* Roe average */
  double roli   = 1.0/rol;
  double rori   = 1.0/ror;
  double rrol   = sqrt(rol);
  double rror   = sqrt(ror);
  double rroi   = 1.0/(rrol+rror);
  double robar  = rrol*rror;
  double robari = 1.0/robar;
  double rrobar = sqrt(robar);
  double rrobai = sqrt(robari);
  double vnbar  = (rrol*vnl+rror*vnr)*rroi;
  double vtbar  = (rrol*vtl+rror*vtr)*rroi;
  double vubar  = (rrol*vul+rror*vur)*rroi;
  double bnbar  = bnc;
  double btbar  = (rrol*btr+rror*btl)*rroi;
  double bubar  = (rrol*bur+rror*bul)*rroi;
  double hl     = (+gamma*gammam1i*prl*roli
		   +0.5*(vnl*vnl+vtl*vtl+vul*vul)
		   +(bnbar*bnbar+btl*btl+bul*bul)*roli);
  double hr     = (+gamma*gammam1i*prr*rori
		   +0.5*(vnr*vnr+vtr*vtr+vur*vur)
		   +(bnbar*bnbar+btr*btr+bur*bur)*rori);
  double hbar   = (rrol*hl+rror*hr)*rroi;
  double btave  = 0.5*(btl+btr);
  double buave  = 0.5*(bul+bur);
  double sgn    = (bnbar > 0)?(1.0):(-1.0);
  /* Characteristic speed */
  double db2 = 0.5*gammam2*gammam1i*((btr-btl)*(btr-btl)+(bur-bul)*(bur-bul))*rroi*rroi;
  double c2  = gammam1*(+hbar
			-0.5*(vnbar*vnbar+vtbar*vtbar+vubar*vubar)-db2
			-(bnbar*bnbar+btbar*btbar+bubar*bubar)*robari);
  double a2  = (gammam1*(+hbar-0.5*(vnbar*vnbar+vtbar*vtbar+vubar*vubar)-db2)
		- gammam2*(bnbar*bnbar+btbar*btbar+bubar*bubar)*robari);
  double ca2 = bnbar*bnbar*robari;
  double cb2 = (btbar*btbar+bubar*bubar)*robari;
  double cf2 = 0.5*(a2+sqrt(cb2*(a2+c2+ca2)+(c2-ca2)*(c2-ca2)));
  double cs2 = c2*ca2/cf2;
  double cf  = sqrt(cf2);
  double cs  = sqrt(cs2);
  double ca  = sqrt(ca2);
  double c   = sqrt(c2);
  double cfi = 1.0/cf;
  /* Renormalization */
  double sgr1 = btbar*btbar+bubar*bubar-FL_EPS;
  double sf11 = (sgr1 > 0)?1.0:0.0;
  double sf21 = 1.0-sf11;
  double beti = 1.0/sqrt(btbar*btbar+bubar*bubar+sf21);
  double bett = sf11*btbar*beti+sf21*sqrt(0.5);
  double betu = sf11*bubar*beti+sf21*sqrt(0.5);
  double sgr2 = (btbar*btbar+bubar*bubar)*robari+fabs(ca2-c2)-FL_EPS;
  double sf12 = (sgr2 > 0)?1.0:0.0;
  double sf22 = 1.0-sf12;
  double cfca = max(0.0,cf2-ca2);
  double cfcs = max(0.0,cf2-cs2);
  double cfa  = max(0.0,cf2-c2);
  double alpi = 1.0/sqrt(cfcs+sf22);
  double alpf = sf12*sqrt(cfca)*alpi+sf22;
  double alps = sf12*sqrt(cfa )*alpi;
  /* Eigenvalues & entropy correction */
  double eeps = 0.0;
  /* eeps = (vnr-vnl+abs(vnr-vnl))*0.25; */
  double evpf = max(fabs(vnbar+cf),eeps);
  double evmf = max(fabs(vnbar-cf),eeps);
  double evps = max(fabs(vnbar+cs),eeps);
  double evms = max(fabs(vnbar-cs),eeps);
  double evpa = max(fabs(vnbar+ca),eeps);
  double evma = max(fabs(vnbar-ca),eeps);
  double ev00 = max(fabs(vnbar),eeps);
  /* Amplitudes */
  double dltro = ror-rol;
  double dltmn = robar*(vnr-vnl);
  double dltmt = robar*(vtr-vtl);
  double dltmu = robar*(vur-vul);
  double dltbt = btr-btl;
  double dltbu = bur-bul;
  double t1    = bett*dltbt+betu*dltbu;
  double t2    = (prr-prl+(btave*dltbt+buave*dltbu)+gammam2*(btbar*dltbt+bubar*dltbu))*gammam1i;
  double t3    = betu*dltbt-bett*dltbu;
  double s1    = dltmn;
  double s2    = bett*dltmt+betu*dltmu;
  double s3    = betu*dltmt-bett*dltmu;
  double p11   = alps*cf*rrobai;
  double p12   =-alpf*c2*cfi*rrobai;
  double p21   = alpf*(cf2-c2*gammam2*gammam1i);
  double p22   = alps*(cs2-c2*gammam2*gammam1i);
  double q11   = alpf*cf;
  double q12   = alps*cs;
  double q21   =-alps*ca*sgn;
  double q22   = alpf*c*sgn;
  double detpi = 1.0/(p11*p22-p12*p21);
  double detqi = 1.0/(q11*q22-q12*q21);
  double wvpf  = 0.5*(( p22*t1-p12*t2)*detpi+( q22*s1-q12*s2)*detqi);
  double wvmf  = 0.5*(( p22*t1-p12*t2)*detpi-( q22*s1-q12*s2)*detqi);
  double wvps  = 0.5*((-p21*t1+p11*t2)*detpi+(-q21*s1+q11*s2)*detqi);
  double wvms  = 0.5*((-p21*t1+p11*t2)*detpi-(-q21*s1+q11*s2)*detqi);
  double wvpa  = 0.5*(rrobar*t3-sgn*s3);
  double wvma  = 0.5*(rrobar*t3+sgn*s3);
  double wv00  = dltro-alpf*(wvpf+wvmf)-alps*(wvps+wvms);
  /* Fluxes at left/right faces */
  double frol = rol*vnl;
  double fmnl = rol*vnl*vnl+prl+0.5*(-bnbar*bnbar+btl*btl+bul*bul);
  double fmtl = rol*vnl*vtl-bnbar*btl;
  double fmul = rol*vnl*vul-bnbar*bul;
  double fbtl = vnl*btl-vtl*bnbar;
  double fbul = vnl*bul-vul*bnbar;
  double fenl = rol*vnl*hl-bnbar*(bnbar*vnl+btl*vtl+bul*vul);
  double fror = ror*vnr;
  double fmnr = ror*vnr*vnr+prr+0.5*(-bnbar*bnbar+btr*btr+bur*bur);
  double fmtr = ror*vnr*vtr-bnbar*btr;
  double fmur = ror*vnr*vur-bnbar*bur;
  double fbtr = vnr*btr-vtr*bnbar;
  double fbur = vnr*bur-vur*bnbar;
  double fenr = ror*vnr*hr-bnbar*(bnbar*vnr+btr*vtr+bur*vur);
  /* Eigenvectors */
  double rvpfro = alpf;
  double rvpfmn = alpf*(vnbar+cf);
  double rvpfmt = alpf*vtbar-alps*bett*ca*sgn;
  double rvpfmu = alpf*vubar-alps*betu*ca*sgn;
  double rvpfbt = alps*bett*cf*rrobai;
  double rvpfbu = alps*betu*cf*rrobai;
  double rvpfen = alpf*(0.5*(vnbar*vnbar+vtbar*vtbar+vubar*vubar)+db2+cf*vnbar+cf2*gammam1i+(cf2-c2)*gammam2*gammam1i)- alps*ca*(bett*vtbar+betu*vubar)*sgn;
  double rvmfro = alpf;
  double rvmfmn = alpf*(vnbar-cf);
  double rvmfmt = alpf*vtbar+alps*bett*ca*sgn;
  double rvmfmu = alpf*vubar+alps*betu*ca*sgn;
  double rvmfbt = rvpfbt;
  double rvmfbu = rvpfbu;
  double rvmfen = alpf*(0.5*(vnbar*vnbar+vtbar*vtbar+vubar*vubar)+db2-cf*vnbar+cf2*gammam1i+(cf2-c2)*gammam2*gammam1i)+alps*ca*(bett*vtbar+betu*vubar)*sgn;
  double rvpsro = alps;
  double rvpsmn = alps*(vnbar+cs);
  double rvpsmt = alps*vtbar+c*sgn*alpf*bett;
  double rvpsmu = alps*vubar+c*sgn*alpf*betu;
  double rvpsbt =-rrobai*c2*alpf*bett*cfi;
  double rvpsbu =-rrobai*c2*alpf*betu*cfi;
  double rvpsen = alps*(0.5*(vnbar*vnbar+vtbar*vtbar+vubar*vubar)+db2+cs*vnbar+cs2*gammam1i+(cs2-c2)*gammam2*gammam1i)+alpf*c*(bett*vtbar+betu*vubar)*sgn;
  double rvmsro = alps;
  double rvmsmn = alps*(vnbar-cs);
  double rvmsmt = alps*vtbar-c*sgn*alpf*bett;
  double rvmsmu = alps*vubar-c*sgn*alpf*betu;
  double rvmsbt = rvpsbt;
  double rvmsbu = rvpsbu;
  double rvmsen = alps*(0.5*(vnbar*vnbar+vtbar*vtbar+vubar*vubar)+db2-cs*vnbar+cs2*gammam1i+(cs2-c2)*gammam2*gammam1i)-alpf*c*(bett*vtbar+betu*vubar)*sgn;
  double rvparo = 0.0;
  double rvpamn = 0.0;
  double rvpamt =-sgn*betu;
  double rvpamu = sgn*bett;
  double rvpabt = rrobai*betu;
  double rvpabu =-rrobai*bett;
  double rvpaen =-(betu*vtbar-bett*vubar)*sgn;
  double rvmaro = 0.0;
  double rvmamn = 0.0;
  double rvmamt =-rvpamt;
  double rvmamu =-rvpamu;
  double rvmabt = rvpabt;
  double rvmabu = rvpabu;
  double rvmaen =-rvpaen;
  double rv00ro = 1.0;
  double rv00mn = vnbar;
  double rv00mt = vtbar;
  double rv00mu = vubar;
  double rv00bt = 0.0;
  double rv00bu = 0.0;
  double rv00en = 0.5*(vnbar*vnbar+vtbar*vtbar+vubar*vubar)+db2;
  /* Numerical fluxes */
  *fro = 0.5*(+frol+fror
	      -evpf*wvpf*rvpfro
	      -evmf*wvmf*rvmfro
	      -evps*wvps*rvpsro
	      -evms*wvms*rvmsro
	      -evpa*wvpa*rvparo
	      -evma*wvma*rvmaro
	      -ev00*wv00*rv00ro);
  *fmn = 0.5*(+fmnl+fmnr
	      -evpf*wvpf*rvpfmn
	      -evmf*wvmf*rvmfmn
	      -evps*wvps*rvpsmn
	      -evms*wvms*rvmsmn
	      -evpa*wvpa*rvpamn
	      -evma*wvma*rvmamn
	      -ev00*wv00*rv00mn);
  *fmt = 0.5*(fmtl+fmtr
	      -evpf*wvpf*rvpfmt
	      -evmf*wvmf*rvmfmt
	      -evps*wvps*rvpsmt
	      -evms*wvms*rvmsmt
	      -evpa*wvpa*rvpamt
	      -evma*wvma*rvmamt
	      -ev00*wv00*rv00mt);
  *fmu = 0.5*(fmul+fmur
	      -evpf*wvpf*rvpfmu
	      -evmf*wvmf*rvmfmu
	      -evps*wvps*rvpsmu
	      -evms*wvms*rvmsmu
	      -evpa*wvpa*rvpamu
	      -evma*wvma*rvmamu
	      -ev00*wv00*rv00mu);
  *fbt = 0.5*(fbtl+fbtr
	      -evpf*wvpf*rvpfbt
	      -evmf*wvmf*rvmfbt
	      -evps*wvps*rvpsbt
	      -evms*wvms*rvmsbt
	      -evpa*wvpa*rvpabt
	      -evma*wvma*rvmabt
	      -ev00*wv00*rv00bt);
  *fbu = 0.5*(fbul+fbur
	      -evpf*wvpf*rvpfbu
	      -evmf*wvmf*rvmfbu
	      -evps*wvps*rvpsbu
	      -evms*wvms*rvmsbu
	      -evpa*wvpa*rvpabu
	      -evma*wvma*rvmabu
	      -ev00*wv00*rv00bu);
  *fen = 0.5*(fenl+fenr
	      -evpf*wvpf*rvpfen
	      -evmf*wvmf*rvmfen
	      -evps*wvps*rvpsen
	      -evms*wvms*rvmsen
	      -evpa*wvpa*rvpaen
	      -evma*wvma*rvmaen
	      -ev00*wv00*rv00en);
}

void calc_flux_hlld(double rol, double vnl, double vtl, double vul, double btl, double bul, double prl,
		    double ror, double vnr, double vtr, double vur, double btr, double bur, double prr,
		    double bnc, double gamma, const double *dv,
		    double *fro, double *fmn, double *fmt, double *fmu, double *fbt, double *fbu, double *fen)
/* Calculate HLLD fluxes */
/* rol,vnl,vtl,vul,btl,bul,prl: input primitive variables at the left side */
/* ror,vnr,vtr,vur,btr,bur,prr: input primitive variables at the right side */
/* bnc: input normal magnetic field at the interface */
/* gamma: specific heat ratio */
/* fro,fmn,fmt,fmu,fbt,fbu,fen: output HLLD fluxes at the interface*/

/* dv is dummy (not used) */
{
  double gammam1i=1.0/(gamma-1.0);
  /* Bn at the interface */
  double bnc2=bnc*bnc;
  double sgn=(bnc > 0)?(1.0):(-1.0);
  /* Variables at the left-face */
  double roli=1.0/rol;
  double pml=0.5*(btl*btl+bul*bul);
  double ptl=prl+pml;
  double enl=gammam1i*prl+pml+0.5*rol*(vnl*vnl+vtl*vtl+vul*vul);
  double vbl=vtl*btl+vul*bul;
  /* Variables at the right-face */
  double rori=1.0/ror;
  double pmr=0.5*(btr*btr+bur*bur);
  double ptr=prr+pmr;
  double enr=gammam1i*prr+pmr+0.5*ror*(vnr*vnr+vtr*vtr+vur*vur);
  double vbr=vtr*btr+vur*bur;
  /* Maximum/minimum wave speeds */
  double cl2=gamma*prl*roli;
  double cr2=gamma*prr*rori;
  double cal2=bnc2*roli;
  double car2=bnc2*rori;
  double cbl2=cl2+cal2+2.0*pml*roli;
  double cbr2=cr2+car2+2.0*pmr*rori;
  double cfl2=0.5*(cbl2+sqrt(fabs(cbl2*cbl2-4.0*cl2*cal2)));
  double cfr2=0.5*(cbr2+sqrt(fabs(cbr2*cbr2-4.0*cr2*car2)));
  double cfl=sqrt(cfl2);
  double cfr=sqrt(cfr2);
  double cmax=max(cfl,cfr);
  double sl=min(0.0,min(vnl,vnr)-cmax);
  double sr=max(0.0,max(vnl,vnr)+cmax);
  /* HLL average of the normal velocity and the total pressure */
  double slvl=sl-vnl;
  double srvr=sr-vnr;
  double rslvl=rol*slvl;
  double rsrvr=ror*srvr;
  double drsvi=1.0/(rsrvr-rslvl);
  double vnc=(rsrvr*vnr-rslvl*vnl-ptr+ptl)*drsvi;
  double ptc=(rsrvr*ptl-rslvl*ptr+rsrvr*rslvl*(vnr-vnl))*drsvi;
  /* Variables of the outer sides in the Riemann fan */
  double slvc=sl-vnc;
  double srvc=sr-vnc;
  double ro2l=rslvl/slvc;
  double ro2r=rsrvr/srvc;
  double rhdl=rslvl*slvc-bnc2;
  double rhdr=rsrvr*srvc-bnc2;
  double vt2l,vu2l,bt2l,bu2l;
  double vt2r,vu2r,bt2r,bu2r;
  if (fabs(rhdl) > FL_EPS){
    double rhdli=1.0/rhdl;
    double rhnvl=(vnl-vnc)*bnc;
    double rhnbl=rslvl*slvl-bnc2;
    vt2l=vtl+rhnvl*rhdli*btl;
    vu2l=vul+rhnvl*rhdli*bul;
    bt2l=rhnbl*rhdli*btl;
    bu2l=rhnbl*rhdli*bul;
  } else{
    vt2l=vtl;
    vu2l=vul;
    bt2l=btl;
    bu2l=bul;
  }
  if (fabs(rhdr) > FL_EPS){
    double rhdri=1.0/rhdr;
    double rhnvr=(vnr-vnc)*bnc;
    double rhnbr=rsrvr*srvr-bnc2;
    vt2r=vtr+rhnvr*rhdri*btr;
    vu2r=vur+rhnvr*rhdri*bur;
    bt2r=rhnbr*rhdri*btr;
    bu2r=rhnbr*rhdri*bur;
  } else{
    vt2r=vtr;
    vu2r=vur;
    bt2r=btr;
    bu2r=bur;
  }
  double vb2l=vt2l*bt2l+vu2l*bu2l;
  double vb2r=vt2r*bt2r+vu2r*bu2r;
  double en2l=(slvl*enl-ptl*vnl+ptc*vnc+bnc*(vbl-vb2l))/slvc;
  double en2r=(srvr*enr-ptr*vnr+ptc*vnc+bnc*(vbr-vb2r))/srvc;
  /* Variables of the inner sides in the Riemann fan */
  double rro2l=sqrt(ro2l);
  double rro2r=sqrt(ro2r);
  double rro2i=1.0/(rro2r+rro2l);
  double vt3m=(rro2r*vt2r+rro2l*vt2l+(bt2r-bt2l)*sgn)*rro2i;
  double vu3m=(rro2r*vu2r+rro2l*vu2l+(bu2r-bu2l)*sgn)*rro2i;
  double bt3m=(rro2l*bt2r+rro2r*bt2l+rro2r*rro2l*(vt2r-vt2l)*sgn)*rro2i;
  double bu3m=(rro2l*bu2r+rro2r*bu2l+rro2r*rro2l*(vu2r-vu2l)*sgn)*rro2i;
  double vb3m=vt3m*bt3m+vu3m*bu3m;
  double en3l=en2l-rro2l*(vb2l-vb3m)*sgn;
  double en3r=en2r+rro2r*(vb2r-vb3m)*sgn;
  /* Variables at the interface */
  double rou,vtu,vuu,btu,buu,enu;
  double hl,hr,h2l,h3l,h2r,h3r;
  hl=(vnc > 0)?(1.0):(0.0);
  hr=1.0-hl;
  h2l=(vnc-fabs(bnc)/rro2l > 0)?(1.0):(0.0);
  h3l=(1.0-h2l)*hl;
  h2r=(vnc+fabs(bnc)/rro2r > 0)?(0.0):(1.0);
  h3r=(1.0-h2r)*hr;
  rou=ro2l*hl+ro2r*hr;
  vtu=vt2l*h2l+vt3m*(h3l+h3r)+vt2r*h2r;
  vuu=vu2l*h2l+vu3m*(h3l+h3r)+vu2r*h2r;
  btu=bt2l*h2l+bt3m*(h3l+h3r)+bt2r*h2r;
  buu=bu2l*h2l+bu3m*(h3l+h3r)+bu2r*h2r;
  enu=en2l*h2l+en3l*h3l+en3r*h3r+en2r*h2r;
  /* HLLD fluxes */
  *fro=rou*vnc;
  *fmn=rou*vnc*vnc+ptc-bnc2*0.5;
  *fmt=rou*vtu*vnc-bnc*btu;
  *fmu=rou*vuu*vnc-bnc*buu;
  *fbt=btu*vnc-bnc*vtu;
  *fbu=buu*vnc-bnc*vuu;
  *fen=(enu+ptc)*vnc-bnc*(vtu*btu+vuu*buu);
}

void calc_flux_lhlld(double rol, double vnl, double vtl, double vul, double btl, double bul, double prl,
		     double ror, double vnr, double vtr, double vur, double btr, double bur, double prr,
		     double bnc, double gamma, const double *dv,
		     double *fro, double *fmn, double *fmt, double *fmu, double *fbt, double *fbu, double *fen)
/* Calculate low-dissipation HLLD fluxes (quasi-all-speed extension) */
/* rol,vnl,vtl,vul,btl,bul,prl: input primitive variables at the left side */
/* ror,vnr,vtr,vur,btr,bur,prr: input primitive variables at the right side */
/* bnc: input normal magnetic field at the interface */
/* gamma: specific heat ratio */
/* fro,fmn,fmt,fmu,fbt,fbu,fen: output HLLD fluxes at the interface*/

/* Shock detection using dv */
{
  double gammam1i=1.0/(gamma-1.0);
  /* Bn at the interface */
  double bnc2=bnc*bnc;
  double sgn=(bnc > 0)?(1.0):(-1.0);
  /* Variables at the left-face */
  double roli=1.0/rol;
  double vl2=vnl*vnl+vtl*vtl+vul*vul;
  double pml=0.5*(btl*btl+bul*bul);
  double ptl=prl+pml;
  double enl=gammam1i*prl+pml+0.5*rol*vl2;
  double vbl=vtl*btl+vul*bul;
  /* Variables at the right-face */
  double rori=1.0/ror;
  double vr2=vnr*vnr+vtr*vtr+vur*vur;
  double pmr=0.5*(btr*btr+bur*bur);
  double ptr=prr+pmr;
  double enr=gammam1i*prr+pmr+0.5*ror*vr2;
  double vbr=vtr*btr+vur*bur;
  /* Maximum/minimum wave speeds */
  double cl2=gamma*prl*roli;
  double cr2=gamma*prr*rori;
  double cal2=bnc2*roli;
  double car2=bnc2*rori;
  double cbl2=cl2+cal2+2.0*pml*roli;
  double cbr2=cr2+car2+2.0*pmr*rori;
  double cfl2=0.5*(cbl2+sqrt(fabs(cbl2*cbl2-4.0*cl2*cal2)));
  double cfr2=0.5*(cbr2+sqrt(fabs(cbr2*cbr2-4.0*cr2*car2)));
  cbl2=vl2+cal2+2.0*pml*roli; /* Sound vel => convective vel. */
  cbr2=vr2+car2+2.0*pmr*rori; /* Sound vel => convective vel. */
  double ccl2=0.5*(cbl2+sqrt(fabs(cbl2*cbl2-4.0*vl2*cal2))); /* Sound vel => convective vel. */
  double ccr2=0.5*(cbr2+sqrt(fabs(cbr2*cbr2-4.0*vr2*car2))); /* Sound vel => convective vel. */
  double cmax=sqrt(max(cfl2,cfr2));
  double sl=min(0.0,min(vnl,vnr)-cmax);
  double sr=max(0.0,max(vnl,vnr)+cmax);
  /* HLL average of the normal velocity and the total pressure */
  double slvl=sl-vnl;
  double srvr=sr-vnr;
  double rslvl=rol*slvl;
  double rsrvr=ror*srvr;
  double drsvi=1.0/(rsrvr-rslvl);
  double theta=min(1.0,(cmax-min(dv[0],0.0))/(cmax-min(dv[1],0.0))); /* Shock detection */
  theta*=theta;
  theta*=theta;
  double vnc=(rsrvr*vnr-rslvl*vnl-theta*(ptr-ptl))*drsvi;
  double mcc=min(1.0,sqrt(max(ccl2,ccr2))/cmax);
  double ptc=(rsrvr*ptl-rslvl*ptr
	      +(mcc*(2.0-mcc))*rsrvr*rslvl*(vnr-vnl))*drsvi; /* Quasi-all-speed extension */
  /* Variables of the outer sides in the Riemann fan */
  double slvc=sl-vnc;
  double srvc=sr-vnc;
  double ro2l=rslvl/slvc;
  double ro2r=rsrvr/srvc;
  double rhdl=rslvl*slvc-bnc2;
  double rhdr=rsrvr*srvc-bnc2;
  double vt2l,vu2l,bt2l,bu2l;
  double vt2r,vu2r,bt2r,bu2r;
  if (fabs(rhdl) > FL_EPS){
    double rhdli=1.0/rhdl;
    double rhnvl=(vnl-vnc)*bnc;
    double rhnbl=rslvl*slvl-bnc2;
    vt2l=vtl+rhnvl*rhdli*btl;
    vu2l=vul+rhnvl*rhdli*bul;
    bt2l=rhnbl*rhdli*btl;
    bu2l=rhnbl*rhdli*bul;
  } else{
    vt2l=vtl;
    vu2l=vul;
    bt2l=btl;
    bu2l=bul;
  }
  if (fabs(rhdr) > FL_EPS){
    double rhdri=1.0/rhdr;
    double rhnvr=(vnr-vnc)*bnc;
    double rhnbr=rsrvr*srvr-bnc2;
    vt2r=vtr+rhnvr*rhdri*btr;
    vu2r=vur+rhnvr*rhdri*bur;
    bt2r=rhnbr*rhdri*btr;
    bu2r=rhnbr*rhdri*bur;
  } else{
    vt2r=vtr;
    vu2r=vur;
    bt2r=btr;
    bu2r=bur;
  }
  double vb2l=vt2l*bt2l+vu2l*bu2l;
  double vb2r=vt2r*bt2r+vu2r*bu2r;
  double en2l=(slvl*enl-ptl*vnl+ptc*vnc+bnc*(vbl-vb2l))/slvc;
  double en2r=(srvr*enr-ptr*vnr+ptc*vnc+bnc*(vbr-vb2r))/srvc;
  /* Variables of the inner sides in the Riemann fan */
  double rro2l=sqrt(ro2l);
  double rro2r=sqrt(ro2r);
  double rro2i=1.0/(rro2r+rro2l);
  double vt3m=(rro2r*vt2r+rro2l*vt2l+(bt2r-bt2l)*sgn)*rro2i;
  double vu3m=(rro2r*vu2r+rro2l*vu2l+(bu2r-bu2l)*sgn)*rro2i;
  double bt3m=(rro2l*bt2r+rro2r*bt2l+rro2r*rro2l*(vt2r-vt2l)*sgn)*rro2i;
  double bu3m=(rro2l*bu2r+rro2r*bu2l+rro2r*rro2l*(vu2r-vu2l)*sgn)*rro2i;
  double vb3m=vt3m*bt3m+vu3m*bu3m;
  double en3l=en2l-rro2l*(vb2l-vb3m)*sgn;
  double en3r=en2r+rro2r*(vb2r-vb3m)*sgn;
  /* Variables at the interface */
  double rou,vtu,vuu,btu,buu,enu;
  double hl,hr,h2l,h3l,h2r,h3r;
  hl=(vnc > 0)?(1.0):(0.0);
  hr=1.0-hl;
  h2l=(vnc-fabs(bnc)/rro2l > 0)?(1.0):(0.0);
  h3l=(1.0-h2l)*hl;
  h2r=(vnc+fabs(bnc)/rro2r > 0)?(0.0):(1.0);
  h3r=(1.0-h2r)*hr;
  rou=ro2l*hl+ro2r*hr;
  vtu=vt2l*h2l+vt3m*(h3l+h3r)+vt2r*h2r;
  vuu=vu2l*h2l+vu3m*(h3l+h3r)+vu2r*h2r;
  btu=bt2l*h2l+bt3m*(h3l+h3r)+bt2r*h2r;
  buu=bu2l*h2l+bu3m*(h3l+h3r)+bu2r*h2r;
  enu=en2l*h2l+en3l*h3l+en3r*h3r+en2r*h2r;
  /* HLLD fluxes */
  *fro=rou*vnc;
  *fmn=rou*vnc*vnc+ptc-bnc2*0.5;
  *fmt=rou*vtu*vnc-bnc*btu;
  *fmu=rou*vuu*vnc-bnc*buu;
  *fbt=btu*vnc-bnc*vtu;
  *fbu=buu*vnc-bnc*vuu;
  *fen=(enu+ptc)*vnc-bnc*(vtu*btu+vuu*buu);
}

void calc_flux_viscous(const double *vx, const double *vy, const double *vz,
		       double rnu, double idx, double idy, double idz,
		       int xoffset, int yoffset, int zoffset,
		       double (*func_fc)(const double *f),
		       double (*func_df)(const double *f),
		       double *fmx, double *fmy, double *fmz, double *fen)
/* Calculate viscous fluxes at the interface i-1/2 */
/* v{x,y,z} are defined @ cell center (i,j,k) */
/* rnu = viscosity coefficient (ro*nu) at the interface, idx=1.0/dx */
/* xoffset=1, yoffset=nx (or 0 in 1D), zoffset=nx*ny (or 0 in 2D) */
/* fmx,fmy,fmz,fen: output fluxes at the interface for momentum and energy */
{
  int ss[4]={-2*xoffset,-xoffset,0,+xoffset}; /* i-2,i-1,i,i+1 */
  double c2=2./3.,c4=4./3.;

  double vxi[4]={vx[ss[0]],vx[ss[1]],vx[ss[2]],vx[ss[3]]};
  double vyi[4]={vy[ss[0]],vy[ss[1]],vy[ss[2]],vy[ss[3]]};
  double vzi[4]={vz[ss[0]],vz[ss[1]],vz[ss[2]],vz[ss[3]]};
  double vx0=func_fc(&vxi[1]); /* vx @ i-1/2 */
  double vy0=func_fc(&vyi[1]); /* vy @ i-1/2 */
  double vz0=func_fc(&vzi[1]); /* vz @ i-1/2 */
  double dvxdx=func_df(&vxi[1])*idx; /* dvx/dx @ i-1/2 */
  double dvydx=func_df(&vyi[1])*idx; /* dvy/dx @ i-1/2 */
  double dvzdx=func_df(&vzi[1])*idx; /* dvz/dx @ i-1/2 */

  double dv[4];
  /* 2nd order difference in orthogonal direction for simplicity */
  dv[0]=0.5*(vx[ss[0]+yoffset]-vx[ss[0]-yoffset]);
  dv[1]=0.5*(vx[ss[1]+yoffset]-vx[ss[1]-yoffset]);
  dv[2]=0.5*(vx[ss[2]+yoffset]-vx[ss[2]-yoffset]);
  dv[3]=0.5*(vx[ss[3]+yoffset]-vx[ss[3]-yoffset]);
  double dvxdy=func_fc(&dv[1])*idy; /* dvx/dy @ i-1/2 */
  dv[0]=0.5*(vy[ss[0]+yoffset]-vy[ss[0]-yoffset]);
  dv[1]=0.5*(vy[ss[1]+yoffset]-vy[ss[1]-yoffset]);
  dv[2]=0.5*(vy[ss[2]+yoffset]-vy[ss[2]-yoffset]);
  dv[3]=0.5*(vy[ss[3]+yoffset]-vy[ss[3]-yoffset]);
  double dvydy=func_fc(&dv[1])*idy; /* dvy/dy @ i-1/2 */
  dv[0]=0.5*(vx[ss[0]+zoffset]-vx[ss[0]-zoffset]);
  dv[1]=0.5*(vx[ss[1]+zoffset]-vx[ss[1]-zoffset]);
  dv[2]=0.5*(vx[ss[2]+zoffset]-vx[ss[2]-zoffset]);
  dv[3]=0.5*(vx[ss[3]+zoffset]-vx[ss[3]-zoffset]);
  double dvxdz=func_fc(&dv[1])*idz; /* dvx/dz @ i-1/2 */
  dv[0]=0.5*(vz[ss[0]+zoffset]-vz[ss[0]-zoffset]);
  dv[1]=0.5*(vz[ss[1]+zoffset]-vz[ss[1]-zoffset]);
  dv[2]=0.5*(vz[ss[2]+zoffset]-vz[ss[2]-zoffset]);
  dv[3]=0.5*(vz[ss[3]+zoffset]-vz[ss[3]-zoffset]);
  double dvzdz=func_fc(&dv[1])*idz; /* dvz/dz @ i-1/2 */

  *fmx=-rnu*(c4*dvxdx-c2*(dvydy+dvzdz));
  *fmy=-rnu*(dvydx+dvxdy);
  *fmz=-rnu*(dvzdx+dvxdz);
  *fen=(*fmx)*vx0+(*fmy)*vy0+(*fmz)*vz0;
}

void calc_flux_resistive(const double *by, const double *bz, const double *ey, const double *ez,
			 int xoffset, int yoffset, int zoffset,
			 double (*func_fc)(const double *f),
			 double *fby, double *fbz, double *fen)
/* Calculate resistive flux at the interface i-1/2 */
/* CT grid spacing is considered. */
/* B are defined @ By(i,j-1/2,k), Bz(i,j,k-1/2) */
/* E are defined @ Ey(i-1/2,j,k-1/2), Ez(i-1/2,j-1/2,k), prepared through MHD_CT_ERES */
/* xoffset=1, yoffset=nx (or 0 in 1D), zoffset=nx*ny (or 0 in 2D) */
/* fby,fbz,fen: output fluxes at the interface for B and energy */
{
  int ss[4]={-2*xoffset,-xoffset,0,+xoffset}; /* i-2,i-1,i,i+1 */
  /* 2nd order average in orthogonal direction for simplicity */
  double ey0=0.5*(ey[0]+ey[zoffset]); /* i-1/2,j,k */
  double ez0=0.5*(ez[0]+ez[yoffset]); /* i-1/2,j,k */
  double byi[4]={0.5*(by[ss[0]]+by[ss[0]+yoffset]),
		 0.5*(by[ss[1]]+by[ss[1]+yoffset]),
		 0.5*(by[ss[2]]+by[ss[2]+yoffset]),
		 0.5*(by[ss[3]]+by[ss[3]+yoffset])}; /* (i-2..i+1),j,k */
  double bzi[4]={0.5*(bz[ss[0]]+bz[ss[0]+zoffset]),
		 0.5*(bz[ss[1]]+bz[ss[1]+zoffset]),
		 0.5*(bz[ss[2]]+bz[ss[2]+zoffset]),
		 0.5*(bz[ss[3]]+bz[ss[3]+zoffset])}; /* (i-2..i+1),j,k */
  double by0=func_fc(&byi[1]);			     /* i-1/2,j,k */
  double bz0=func_fc(&bzi[1]);			     /* i-1/2,j,k */
  *fby=-ez0;					     /* Not used for CT-By */
  *fbz=+ey0;					     /* Not used for CT-Bz */
  *fen=ey0*bz0-ez0*by0;
}

void calc_flux_hall_lf(double rol, double vnl, double vtl, double vul, double btl, double bul, double enl,
		       double ror, double vnr, double vtr, double vur, double btr, double bur, double enr,
		       double bnc, double smax,
		       double *fbt, double *fbu, double *fen)
/* Calculate Lax-Friedrichs fluxes for Hall term */
/* rol,vnl,vtl,vul,btl,bul,enl: input variables at the left side. v is Hall velocity -j/nec*/
/* ror,vnr,vtr,vur,btr,bur,enr: input variables at the right side */
/* bnc: input normal magnetic field at the interface */
/* smax: Maximum velocity */
/* fbt,fbu,fen: output LF fluxes at the interface*/
{
  /* Variables at the left-face */
  double pml=0.5*(bnc*bnc+btl*btl+bul*bul);
  double vbl=vnl*bnc+vtl*btl+vul*bul;
  /* Variables at the right-face */
  double pmr=0.5*(bnc*bnc+btr*btr+bur*bur);
  double vbr=vnr*bnc+vtr*btr+vur*bur;
  /* Fluxes at the left-face */
  double fbtl=btl*vnl-bnc*vtl;
  double fbul=bul*vnl-bnc*vul;
  double fenl=2.0*pml*vnl-bnc*vbl;
  /* FLuxes at the right-face */
  double fbtr=btr*vnr-bnc*vtr;
  double fbur=bur*vnr-bnc*vur;
  double fenr=2.0*pmr*vnr-bnc*vbr;
  /* LF fluxes */
  *fbt=0.5*(fbtl+fbtr-smax*(btr-btl));
  *fbu=0.5*(fbul+fbur-smax*(bur-bul));
  *fen=0.5*(fenl+fenr-smax*(pmr-pml));
}
