;+
; NAME:
;   hs_biasgen
;
; PURPOSE:
;   Create a super bias image from all of the single bias images
;
; CALLING SEQUENCE:
;  hs_biasgen, biasfiles, outdir=outdir, mjd=mjd
;
; INPUTS:
;   biasfiles   - a linst of bias files that are to be comined
;
; OPTIONAL KEYWORDS:
;   outdir     - Output directory for the final combined arc
;   mjd        - Integer part of the MJD for the observations
;                this will create the outname for the file if
;                 in the headers of the file
;
; OUTPUTS:
;      creates the /calibration/XXXX/bias-YYYY.fits file where
;      XXXX is the rerun number and YYYY is the MJD
;
; OPTIONAL OUTPUTS:
;
;
; COMMENTS:
;
;
; EXAMPLES:
;
; BUGS:
;
;    Not sure the scaling for the images is appropriate
;    is the entrie chip being used or just the illuminated
;    region
;;
; REVISION HISTORY:
;   April 2004 - Written by R Cool - UofA
;                adapted from code in Schlegals IDLSPEC2d
;                 package
;-
;------------------------------------------------------------------------------



;;This program will combine the biases for the data.

pro hs_biasgen1, files, outfile=outfile, outdir=outdir
  
  nfile = n_elements(files)
  if (NOT keyword_set(sigrej)) then begin
     if (nfile) LE 2 then sigrej = 1.0 $
     else if (nfile EQ 3) then sigrej=1.1 $
     else if (nfile EQ 4) then sigrej = 1.3 $
     else if (nfile EQ 5) then sigrej = 1.6 $
     else if (nfile EQ 6) then sigrej = 1.9 $
     else sigrej = 2.0
  endif
  
  if (NOT keyword_set(maxiter)) then maxiter = 10.0

  for iamp = 1, 4 do begin
  for ifile =0, nfile-1 do begin
     
        splog, 'Reading file #', ifile+1, ' of ', nfile 
        thisimage = readfits(files[ifile], exten_no=iamp, header, /silent)
        
            ;;Overscan and trim correct the bias image
        hs_overscan, thisimage, header
        hs_trim, thisimage, header=header
        thisivar = 1/thisimage
                    
        ;;For the first image, setup the arrays for the combination 
        if (ifile EQ 0) then begin
           hdr0 = header
           imgarr = make_array(dimension=[size(thisimage, $
                                /dimens), nfile], /float)
           inmask = make_array(dimension=$
                               [size(thisimage, /dimens), nfile], /byte)
           sxaddpar, hdr0, 'NEDP', nfile, 'Number of exposures in this file'
        endif
        
        imgarr[*,*,ifile]= thisimage
        inmask[*,*,ifile] = thisivar LT 0
        
     endfor
     
     mnimg = djs_avsigclip(imgarr, sigrej=sigrej, maxiter=maxiter, $
                           inmask=inmask)
     splog, 'Writing combined bias ', $
        djs_filepath(outfile, root_dir = outdir), $
        ' amplifier', string(iamp)
     
     sxdelpar, hdr0, 'EXTNAME'
     sxdelpar, hdr0, 'BSCALE'
     sxdelpar, hdr0, 'BZERO'
     
     if (iamp eq 1) then begin
        
        junk  = headfits(files(0), exten=0)
        sxaddpar, junk, 'bitpix', long(-32)
        writefits, djs_filepath(outfile, $
                                root_dir = outdir), 1, junk
        writefits, djs_filepath(outfile, $
                                root_dir = outdir), mnimg, hdr0, /append
     endif else begin
        writefits, djs_filepath(outfile, $
                                root_dir = outdir), mnimg, hdr0, /append
     endelse 
     
  endfor
  
   
  return
  
end

  

PRO hs_biasgen, biasfiles,  outdir=outdir, indir=indir, mjd=mjd
  
  
  if not keyword_set(outdir) then outdir = cwd()
  if not keyword_set(indir) then indir = cwd()
  
  ;;Grab what we need out of the FITS headers
  
  
  ;;Setup the subregions of the chips
  chip = findgen(4)+1
  stop = 'go'
  imagetyp = strarr(n_elements(biasfiles))
  
  
  
  ;;Check to make sure that the files are all biases
  for ifile = 0, n_elements(biasfiles) - 1l do begin
     
     head = headfits(biasfiles(ifile), exten=1)

     imtype = strmid((sxpar(head, 'imagetyp')),0, 4)     
  endfor
  

  ;;Check to be sure that there are proper numbers of biases in the directory
  if n_elements(biasfiles) GT 50 then begin
     splog, 'Too many biases (' + strn(n_elements(biasfiles)) + $
        ') to load into memory'
     stop = 'stop'
  endif
  
  
  if n_elements(biasfiles) LT 3 then begin
     splog, 'Expected at least 3 biases, got ' + $ 
        strn(n_elements(biasfiles))
     stop = 'stop'
  endif
  
  
  if stop ne 'stop' then begin
     biasname = 'bias.fits'
     splog, 'Generating Bias ' + biasname
     splog, 'Output directory ' + outdir
     hs_biasgen1, biasfiles, outfile=biasname, outdir=outdir
     
  endif
  
  return
  
end
;--------------------------------------------------------------
     
     
