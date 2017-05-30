PRO fit

  ;RESOLVE_ROUTINE, '/Users/Home/Documents/Work/idl/InPath/pahfit/pahfit'

  plot, [1,1], [1,1], $
        xra=[1,25], /xst, $
        yra=[0,1.1], /yst, /nodata
  
  files = file_search('SB*')
  for i=0, n_elements(files)-1 DO BEGIN
     readcol, files[i], wav, flux
     o = WHERE(wav GE 6. AND wav LT 35.)
     fit=pahfit(wav[o], flux[o], flux[o]/10)

     mn=1.
     mx=25.
     lam=mn+findgen(2400)/(2400)*(mx-mn)

     yfit=pahfit_components(lam,fit,DUST_CONTINUUM=dusts, $
                            TOTAL_DUST_CONTINUUM=dust_tot,STARLIGHT=stars, $
                            DUST_FEATURES=features, $
                            TOTAL_DUST_FEATURES=features_tot, $
                            LINES=lines,TOTAL_LINES=lines_tot, $
                            EXTINCTION_FAC=ext,_EXTRA=e)
     rat = fltarr(n_elements(lam))
     for j=0,(size(features,/DIMENSIONS))[1]-1 do begin 
        rat+=(features[*,j])/lam
     endfor 

  oplot, lam, rat/max(rat)
     if n_elements(allrat) EQ 0 then allrat = fltarr(n_elements(lam), 5)
     allrat[*,i] = rat/max(rat)
  
  endfor

  resistant_mean, allrat, 2., mean, dim=2
  oplot, lam, mean              ;, color=fsc_color('Red')

  forprint, lam, mean, textout='avg_pah.dat'
  
END
