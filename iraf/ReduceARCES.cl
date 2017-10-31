setinst echelle
hedit *.fits dispaxis 1 add=yes verify=no show=yes update=yes
!ls *fits | grep -v 'bias' > inlist
!sed s/.fits/.c.fits/ inlist > outlist
cosmicrays (input="@inlist", output="@outlist", answer="yes", crmasks="", threshold=50.0, fluxratio=10.0, npasses=20, window="5", interactive=no, train=no,objects="", savefile="", plotfile="")
!ls bias*fits > biaslist
zerocombine (input="@biaslist", output="Zero.fits", combine="average", reject="avsigclip", ccdtype="", process=no, delete=no, clobber=no, scale="none", nlow=0, nhigh=1, nkeep=1, mclip=yes, lsigma=3.0, hsigma=3.0, rdnoise="7", gain="3.8", snoise="0.")
!cp outlist inlist
!sed s/.c.fits/.pc.fits/ inlist > outlist
ccdproc (images="@inlist", output="@outlist", ccdtype="", fixpix = yes, overscan=no, trim=yes, zerocor=yes, darkcor=no, flatcor=no, illumcor=no, fringecor=no, readcor=no, scancor=no, readaxis="line", fixfile = "database/badpix.txt", biassec="", trimsec="[200:1850,1:2048]", zero="Zero.fits", dark="", flat="", illum="", fringe="", minreplace=1.0, scantype="shortscan", nscan=1, interactive=no, function="legendre", order=3, naverage=1, niterate=3, low_reject = 3.0, high_reject = 3.0, grow=0.0)
!cp outlist inlist
!sed s/.pc/.cpc/ inlist > outlist
cosmicrays (input="@inlist", output="@outlist", answer="yes", crmasks="", threshold=50.0, fluxratio=10.95, npasses=20, window="5", interactive=no, train=no,objects="", savefile="", plotfile="")
!ls redflat*cpc.fits > redflatlist
!ls blueflat*cpc.fits > blueflatlist
flatcombine (input="@redflatlist", output="redflat", combine="median", reject="avsigclip", ccdtype="", process=no, subsets=no, delete=no, clobber=no, scale="mode", nlow=1, nhigh=1, nkeep=1, mclip=yes, lsigma=3.0, hsigma=3.0, rdnoise="0", gain="1.", snoise="0.", pclip=-0.5)
flatcombine (input="@blueflatlist", output="blueflat", combine="median", reject="avsigclip", ccdtype="", process=no, subsets=no, delete=no, clobber=no, scale="mode", nlow=1, nhigh=1, nkeep=1, mclip=yes, lsigma=3.0, hsigma=3.0, rdnoise="0", gain="1.", snoise="0.", pclip=-0.5)
imarith blueflat + redflat junk
imarith junk / 2.0 superflat
imdel junk
magnify (input = "superflat", output="superflatmag", xmag = 1.0, ymag=4.0, x1=INDEF, x2=INDEF, dx=INDEF, y1=INDEF, y2=INDEF, dy=INDEF, interpolatio="linear", boundary="nearest", constant=0.0, fluxconserve=yes)
hedit (images="superflatmag", fields="CCDSEC", value="[200:1850,1:8189]",add=yes, addonly=yes, delete=no, verify=no, show=yes, update=yes)
apall(input="superflatmag",output="superflatmagr",format="echelle", references="echtrace130522",interactive=no, find=no, recenter=yes, resize=yes, edit=no, trace=yes, fittrace=no, extract=yes, extras=no, review=no, line=INDEF, nsum=10, lower=-14.0, upper=14.0, b_function="chebyshev", b_order=2, b_sample="-22:-15,15:22", b_naverage=-3, b_niterate=3, b_low_reject=3.0, b_high_rejec = 3.0, b_grow = 0.0, width=18.0, radius=18.0, threshold=0.0, minsep=5.0, maxsep=100000.0, order="increasing", aprecenter="", npeaks=INDEF, shift=no, llimit=INDEF, ulimit=INDEF, ylevel=0.05, peak=yes, bkg=yes, r_grow=0.0, avglimits=no, t_nsum=10, t_step=10, t_nlost=3, t_function="legendre", t_order=10, t_sample="1:1651", t_naverage=1, t_niterate=0, t_low_reject=3.0, t_high_reject=3.0, t_grow=0.0, background="fit", skybox=1, weights="variance", pfit="fit1d", clean=yes, saturation=35000.0, readnoise="RDNOISE", gain="GAIN", lsigma=4.0, usigma=4.0, nsubaps=1)
sfit (input="superflatmagr", output="normflat", ask="yes", lines="*", bands="1", type="ratio", replace=no, wavescale=no, logscale=no, override=yes, listonly=no, interactive=no, sample="*", naverage=-5, function="spline3", order=9, low_reject=3.5, high_reject=3.5, niterate=10, grow=1.0)
!grep -v flat outlist > inlist
!sed s/.cpc/.mcpc/ inlist > outlist
magnify (input = "@inlist", output = "@outlist", xmag = 1.0, ymag = 4.0, x1=INDEF, x2=INDEF, dx=INDEF, y1=INDEF, y2=INDEF, dy=INDEF, interpolatio="linear", boundary="nearest", constant=0.0, fluxconserve=yes)
hedit (images="@outlist", fields="CCDSEC", value="[200:1850,1:8189]",add=yes, addonly=yes, delete=no, verify=no, show=yes, update=yes)
!grep ThAr outlist > tharlist
!sed s/mcpc/rmcpc/ tharlist > tharoutlist
apall (input = "@tharlist", output = "@tharoutlist", format = "echelle", references="superflatmag",interactive=no, find=no, recenter=no, resize=no, edit=no, trace=yes, fittrace=no, extract=yes, extras=no, review=no, line=INDEF, nsum=10, lower=-14.0, upper=14.0, b_function="chebyshev", b_order=2, b_sample="-22:-15,15:22", b_naverage=-3, b_niterate=3, b_low_reject=3.0, b_high_rejec = 3.0, b_grow = 0.0, width=18.0, radius=18.0, threshold=0.0, minsep=5.0, maxsep=100000.0, order="increasing", aprecenter="", npeaks=INDEF, shift=no, llimit=INDEF, ulimit=INDEF, ylevel=0.05, peak=yes, bkg=yes, r_grow=0.0, avglimits=no, t_nsum=10, t_step=10, t_nlost=3, t_function="legendre", t_order=10, t_sample="1:1651", t_naverage=1, t_niterate=0, t_low_reject=3.0, t_high_reject=3.0, t_grow=0.0, background="none", skybox=1, weights="none", pfit="fit1d", clean=yes, saturation=35000.0, readnoise="RDNOISE", gain="GAIN", lsigma=4.0, usigma=4.0, nsubaps=1)
!grep -v ThAr outlist > inlist
!sed s/mcpc/rmcpc/ inlist > outlist
apall (input = "@inlist", output = "@outlist", format = "echelle", references="superflatmag",interactive=no, find=no, recenter=no, resize=no, edit=no, trace=no, fittrace=no, extract=no, extras=no, review=no, line=INDEF, nsum=10, lower=-14.0, upper=14.0, b_function="chebyshev", b_order=2, b_sample="-22:-15,15:22", b_naverage=-3, b_niterate=3, b_low_reject=3.0, b_high_rejec = 3.0, b_grow = 0.0, width=18.0, radius=18.0, threshold=0.0, minsep=5.0, maxsep=100000.0, order="increasing", aprecenter="", npeaks=INDEF, shift=no, llimit=INDEF, ulimit=INDEF, ylevel=0.05, peak=yes, bkg=yes, r_grow=0.0, avglimits=no, t_nsum=10, t_step=10, t_nlost=3, t_function="legendre", t_order=10, t_sample="1:1651", t_naverage=1, t_niterate=0, t_low_reject=3.0, t_high_reject=3.0, t_grow=0.0, background="fit", skybox=1, weights="variance", pfit="fit1d", clean=yes, saturation=35000.0, readnoise="RDNOISE", gain="GAIN", lsigma=4.0, usigma=4.0, nsubaps=1)
!sed s/mcpc/noscat/ inlist > noscatlist
imcopy @inlist @noscatlist
imdel @inlist
apscat1 (function = "spline3", order=25, sample="*", naverage=1, low_reject=6.0, high_reject=1.5, niterate=5)
apscat2 (function = "spline3", order=5, sample="*", naverage=1, low_reject=3.0, high_reject=3.0, niterate=2)
apscatter (input="@noscatlist", output="@inlist", apertures="", scatter="", references="@inlist", interactive=no, find=no, recenter=no, resize=no, edit=no, trace=no, fittrace=no, subtract=yes, smooth=yes, fitscatter=yes, fitsmooth=yes, line=INDEF, nsum=-10, buffer=0.4, apscat1="", apscat2='')
apall (input="@inlist", output="@outlist", format = "echelle", references="",profiles="", interactive=no, find=no, recenter=no, resize=no, edit=no, trace=no, fittrace=no, extract=yes, extras=no, review=no, line=INDEF, nsum=10, lower=-14.0, upper=14.0, b_function="chebyshev", b_order=2, b_sample="-22:-15,15:22", b_naverage=-3, b_niterate=3, b_low_reject=3.0, b_high_rejec = 3.0, b_grow = 0.0, width=18.0, radius=18.0, threshold=0.0, minsep=5.0, maxsep=100000.0, order="increasing", aprecenter="", npeaks=INDEF, shift=no, llimit=INDEF, ulimit=INDEF, ylevel=0.05, peak=yes, bkg=yes, r_grow=0.0, avglimits=no, t_nsum=10, t_step=10, t_nlost=3, t_function="legendre", t_order=10, t_sample="1:1651", t_naverage=1, t_niterate=0, t_low_reject=3.0, t_high_reject=3.0, t_grow=0.0, background="fit", skybox=1, weights="variance", pfit="fit1d", clean=yes, saturation=35000.0, readnoise="RDNOISE", gain="GAIN", lsigma=4.0, usigma=4.0, nsubaps=1)
!cp outlist inlist
!sed s/rmcpc/frmcpc/ inlist > outlist
imarith @inlist / normflat @outlist
!mv tharoutlist tharlist
ecreidentify (images="@tharlist", reference="arcnewref.ec", shift=INDEF, cradius=2.0, threshold=50.0, refit=yes, database="database")
!cp outlist inlist
!sed s/frmcpc/wfrmcpc/ inlist > outlist
refspectra (input="@inlist", answer="YES", references="@tharlist", apertures="", refaps="", ignoreaps = yes, select="nearest", sort="date", group="", time=yes, timewrap=17.0, override=yes, confirm=no, assign=yes)
dispcor (input="@inlist", output="@outlist", linearize=no, database="database", w1=INDEF, w2=INDEF, dw=INDEF, nw=INDEF, flux=no, blank=0.0, samedisp=no, global=no, ignoreaps=no, confirm=no, listonly=no, verbose=yes)
setjd (images = "@outlist", observatory="apo", date="date-obs", time="ut", exposure="exptime", ra="ra", dec="dec", epoch="equinox", jd="jd", hjd="hjd", ljd="ljd", utdate=yes, uttime=yes)