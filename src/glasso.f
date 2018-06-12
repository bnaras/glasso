c     mortran 2.0     (version of 7/04/75 mod 7/4/87 (ajc))              
      subroutine glasso(nn,sss,rrho,ia,is,itr,ipen,thr,maxit,www,wwwi,nn    134 
     *iter,ddel,jerr)
      implicit double precision(a-h,o-z)                                    135
      double precision sss(nn,nn),rrho(nn,nn),www(nn,nn),wwwi(nn,nn)        136
      double precision, dimension (:), allocatable :: ss,rho,ww,wwi             
      integer, dimension (:), allocatable :: ir,ie                              
      integer, dimension (:,:), allocatable :: ic                               
      if(ia .eq. 0)goto 10021                                               142
      call lasinv1(nn,sss,rrho,ia,is,itr,ipen,thr,maxit,www,wwwi,nniter,    143 
     *ddel,jerr)
      return                                                                144
10021 continue                                                              145
      allocate(ic(1:2,1:nn),stat=jerr)                                          
      if (jerr.ne.0) return                                                 149
      allocate(ir(1:nn),stat=jerr)                                          150
      if (jerr.ne.0) return                                                 151
      allocate(ie(1:nn),stat=jerr)                                          152
      if(jerr.ne.0) return                                                  153
      call connect(nn,sss,rrho,nc,ic,ir,ie)                                 154
      nnq=0                                                                 154
10030 do 10031 kc=1,nc                                                      154
      nnq=max(ic(2,kc)-ic(1,kc)+1,nnq)                                      154
10031 continue                                                              154
10032 continue                                                              154
      nnq=nnq**2                                                            155
      allocate(ss(1:nnq),stat=jerr)                                         156
      if(jerr.ne.0) return                                                  157
      allocate(rho(1:nnq),stat=jerr)                                        158
      if(jerr.ne.0) return                                                  159
      allocate(ww(1:nnq),stat=jerr)                                         160
      if(jerr.ne.0) return                                                  161
      allocate(wwi(1:nnq),stat=jerr)                                        162
      if(jerr.ne.0) return                                                  163
      nniter=0                                                              163
      ddel=0.0                                                              163
      l=0                                                                   164
10040 do 10041 kc=1,nc                                                      164
      n=ic(2,kc)-ic(1,kc)+1                                                 165
      if(n .gt. 1)goto 10061                                                165
      k=ir(ic(1,kc))                                                        166
      www(:,k)=0.0                                                          166
      www(k,:)=0.0                                                          166
      wwwi(:,k)=0.0                                                         166
      wwwi(k,:)=0.0                                                         166
      goto 10041                                                            167
10061 continue                                                              168
      kb=ic(1,kc)                                                           168
      ke=ic(2,kc)                                                           168
      l=0                                                                   169
10070 do 10071 k=kb,ke                                                      169
      ik=ir(k)                                                              169
10080 do 10081 j=kb,ke                                                      169
      ij=ir(j)                                                              169
      l=l+1                                                                 170
      ss(l)=sss(ij,ik)                                                      170
      rho(l)=rrho(ij,ik)                                                    171
      ww(l)=www(ij,ik)                                                      171
      wwi(l)=wwwi(ij,ik)                                                    172
10081 continue                                                              172
10082 continue                                                              172
10071 continue                                                              173
10072 continue                                                              173
      call lasinv1(n,ss,rho,ia,is,itr,ipen,thr,maxit,ww,wwi,niter,del,je    174 
     *rr)
      if(jerr.ne.0) return                                                  174
      nniter=nniter+niter                                                   174
      ddel=ddel+del                                                         175
10090 do 10091 j=kb,ke                                                      175
      k=ir(j)                                                               176
      www(:,k)=0.0                                                          176
      www(k,:)=0.0                                                          176
      wwwi(:,k)=0.0                                                         176
      wwwi(k,:)=0.0                                                         177
10091 continue                                                              178
10092 continue                                                              178
      l=0                                                                   178
10100 do 10101 k=kb,ke                                                      178
      ik=ir(k)                                                              178
10110 do 10111 j=kb,ke                                                      178
      l=l+1                                                                 178
      wwwi(ir(j),ik)=wwi(l)                                                 178
10111 continue                                                              178
10112 continue                                                              178
10101 continue                                                              179
10102 continue                                                              179
      if(ia .ne. 0)goto 10131                                               179
      l=0                                                                   180
10140 do 10141 k=kb,ke                                                      180
      ik=ir(k)                                                              180
10150 do 10151 j=kb,ke                                                      180
      l=l+1                                                                 180
      www(ir(j),ik)=ww(l)                                                   180
10151 continue                                                              180
10152 continue                                                              180
10141 continue                                                              181
10142 continue                                                              181
10131 continue                                                              182
10041 continue                                                              183
10042 continue                                                              183
      ddel=ddel/nc                                                          184
      if(ia .ne. 0)goto 10171                                               185
10180 do 10181 j=1,nn                                                       185
      if(www(j,j).ne.0.0)goto 10181                                         186
      if(ipen .ne. 0)goto 10201                                             186
      www(j,j)=sss(j,j)                                                     186
      goto 10211                                                            186
10201 continue                                                              186
      www(j,j)=sss(j,j)+rrho(j,j)                                           186
10211 continue                                                              187
10191 continue                                                              187
      wwwi(j,j)=1.0/www(j,j)                                                188
10181 continue                                                              189
10182 continue                                                              189
10171 continue                                                              190
      return                                                                191
      end                                                                   192
      subroutine connect(n,ss,rho,nc,ic,ir,ie)                              193
      implicit double precision(a-h,o-z)                                    194
      double precision ss(n,n),rho(n,n)                                     194
      integer ic(2,n),ir(n),ie(n)                                           195
      ie=0                                                                  195
      nc=0                                                                  195
      is=1                                                                  196
10220 do 10221 k=1,n                                                        196
      if(ie(k).gt.0)goto 10221                                              197
      ir(is)=k                                                              197
      nc=nc+1                                                               197
      ie(k)=nc                                                              197
      ic(1,nc)=is                                                           197
      is=is+1                                                               198
      call row(nc,1,ir((is-1):n),n,ss,rho,ie,na,ir(is:n))                   199
      if(na .ne. 0)goto 10241                                               199
      ic(2,nc)=is-1                                                         199
      goto 10221                                                            199
10241 continue                                                              200
10250 continue                                                              200
10251 continue                                                              200
      nas=na                                                                200
      iss=is                                                                200
      il=iss+nas-1                                                          201
      if(il.ge.n)goto 10252                                                 201
      is=is+na                                                              202
      call row(nc,nas,ir(iss:n),n,ss,rho,ie,na,ir(is:n))                    203
      if(na.eq.0)goto 10252                                                 203
      goto 10251                                                            204
10252 continue                                                              204
      ic(2,nc)=il                                                           205
10221 continue                                                              206
10222 continue                                                              206
      return                                                                207
      end                                                                   208
      subroutine row(nc,nr,jr,n,ss,rho,ie,na,kr)                            209
      implicit double precision(a-h,o-z)                                    210
      double precision ss(n,n),rho(n,n)                                     210
      integer jr(nr),ie(n),kr(*)                                            211
      na=0                                                                  212
10260 do 10261 l=1,nr                                                       212
      k=jr(l)                                                               213
10270 do 10271 j=1,n                                                        213
      if(ie(j).gt.0)goto 10271                                              213
      if(j.eq.k)goto 10271                                                  214
      if(abs(ss(j,k)).le.rho(j,k))goto 10271                                215
      na=na+1                                                               215
      kr(na)=j                                                              215
      ie(j)=nc                                                              216
10271 continue                                                              217
10272 continue                                                              217
10261 continue                                                              218
10262 continue                                                              218
      return                                                                219
      end                                                                   220
      subroutine lasinv1(n,ss,rho,ia,is,itr,ipen,thr,maxit,ww,wwi,niter,    221 
     *del,jerr)
      implicit double precision(a-h,o-z)                                    222
      parameter(eps=1.0e-7)                                                 223
      double precision ss(n,n),rho(n,n),ww(n,n),wwi(n,n)                    224
      double precision, dimension (:,:), allocatable :: vv,xs                   
      double precision, dimension (:), allocatable :: s,x,z,ws,ro,so            
      integer, dimension (:), allocatable :: mm                                 
      nm1=n-1                                                                   
      allocate(vv(1:nm1,1:nm1),stat=jerr)                                       
      if(jerr.ne.0) return                                                      
      if(ia.eq.0) allocate(xs(1:nm1,1:n),stat=jerr)                             
      if(jerr.ne.0) return                                                      
      allocate(s(1:nm1),stat=jerr)                                          235
      if(jerr.ne.0) return                                                  236
      allocate(so(1:nm1),stat=jerr)                                         237
      if(jerr.ne.0) return                                                  238
      allocate(x(1:nm1),stat=jerr)                                          239
      if(jerr.ne.0) return                                                  240
      allocate(z(1:nm1),stat=jerr)                                          241
      if(jerr.ne.0) return                                                  242
      allocate(mm(1:nm1),stat=jerr)                                         243
      if(jerr.ne.0) return                                                  244
      allocate(ro(1:nm1),stat=jerr)                                         245
      if(ia .ne. 0)goto 10291                                               245
      allocate(ws(1:n),stat=jerr)                                           245
10291 continue                                                              246
      if(jerr.ne.0) return                                                  247
      shr=0.0                                                               247
10300 do 10301 j=1,n                                                        247
10310 do 10311 k=1,n                                                        247
      if(j.eq.k)goto 10311                                                  247
      shr=shr+abs(ss(j,k))                                                  247
10311 continue                                                              247
10312 continue                                                              247
10301 continue                                                              248
10302 continue                                                              248
      if(shr .ne. 0.0)goto 10331                                            248
      ww=0.0                                                                248
      wwi=0.0                                                               249
10340 do 10341 j=1,n                                                        249
      if(ipen .ne. 0)goto 10361                                             249
      ww(j,j)=ss(j,j)                                                       249
      goto 10371                                                            249
10361 continue                                                              249
      ww(j,j)=ss(j,j)+rho(j,j)                                              249
10371 continue                                                              250
10351 continue                                                              250
      wwi(j,j)=1.0/max(ww(j,j),eps)                                         251
10341 continue                                                              252
10342 continue                                                              252
      return                                                                253
10331 continue                                                              254
      shr=thr*shr/nm1                                                       255
      if(ia .eq. 0)goto 10391                                               255
      if(is.eq.0) wwi=0.0                                                   256
10400 do 10401 m=1,n                                                        256
      call setup(m,n,ss,rho,ss,vv,s,ro)                                     256
      l=0                                                                   257
10410 do 10411 j=1,n                                                        257
      if(j.eq.m)goto 10411                                                  257
      l=l+1                                                                 257
      x(l)=wwi(j,m)                                                         257
10411 continue                                                              258
10412 continue                                                              258
      call lasso(ro,nm1,vv,s,shr/n,x,z,mm)                                  259
      l=0                                                                   259
10420 do 10421 j=1,n                                                        259
      if(j.eq.m)goto 10421                                                  259
      l=l+1                                                                 259
      wwi(j,m)=x(l)                                                         259
10421 continue                                                              260
10422 continue                                                              260
10401 continue                                                              261
10402 continue                                                              261
      niter=1                                                               261
      return                                                                262
10391 continue                                                              263
      if(is .ne. 0)goto 10441                                               263
      ww=ss                                                                 263
      xs=0.0                                                                263
      goto 10451                                                            264
10441 continue                                                              265
10460 do 10461 j=1,n                                                        265
      xjj=-wwi(j,j)                                                         265
      l=0                                                                   266
10470 do 10471 k=1,n                                                        266
      if(k.eq.j)goto 10471                                                  266
      l=l+1                                                                 266
      xs(l,j)=wwi(k,j)/xjj                                                  266
10471 continue                                                              267
10472 continue                                                              267
10461 continue                                                              268
10462 continue                                                              268
10451 continue                                                              269
10431 continue                                                              269
10480 do 10481 j=1,n                                                        269
      if(ipen .ne. 0)goto 10501                                             269
      ww(j,j)=ss(j,j)                                                       269
      goto 10511                                                            269
10501 continue                                                              269
      ww(j,j)=ss(j,j)+rho(j,j)                                              269
10511 continue                                                              269
10491 continue                                                              269
10481 continue                                                              270
10482 continue                                                              270
      niter=0                                                               271
10520 continue                                                              271
10521 continue                                                              271
      dlx=0.0                                                               272
10530 do 10531 m=1,n                                                        273
      if(itr .eq. 0)goto 10551                                              273
      call intpr('m',1,m,1)                                                 273
10551 continue                                                              274
      x=xs(:,m)                                                             274
      ws=ww(:,m)                                                            275
      call setup(m,n,ss,rho,ww,vv,s,ro)                                     275
      so=s                                                                  276
      call lasso(ro,nm1,vv,s,shr/sum(abs(vv)),x,z,mm)                       277
      l=0                                                                   278
10560 do 10561 j=1,n                                                        278
      if(j.eq.m)goto 10561                                                  279
      l=l+1                                                                 279
      ww(j,m)=so(l)-s(l)                                                    279
      ww(m,j)=ww(j,m)                                                       280
10561 continue                                                              281
10562 continue                                                              281
      dlx=max(dlx,sum(abs(ww(:,m)-ws)))                                     282
      xs(:,m)=x                                                             283
10531 continue                                                              284
10532 continue                                                              284
      niter=niter+1                                                         284
      if(niter.ge.maxit)goto 10522                                          285
      if(dlx.lt.shr)goto 10522                                              285
      goto 10521                                                            286
10522 continue                                                              286
      del=dlx/nm1                                                           286
      call inv(n,ww,xs,wwi)                                                 287
      return                                                                288
      end                                                                   289
      subroutine setup(m,n,ss,rho,ww,vv,s,r)                                290
      implicit double precision(a-h,o-z)                                    291
      double precision ss(n,n),rho(n,n),ww(n,n),vv(n-1,n-1),s(n-1),r(n-1    292 
     *)
      l=0                                                                   293
10570 do 10571 j=1,n                                                        293
      if(j.eq.m)goto 10571                                                  293
      l=l+1                                                                 293
      r(l)=rho(j,m)                                                         293
      s(l)=ss(j,m)                                                          293
      i=0                                                                   294
10580 do 10581 k=1,n                                                        294
      if(k.eq.m)goto 10581                                                  294
      i=i+1                                                                 294
      vv(i,l)=ww(k,j)                                                       294
10581 continue                                                              295
10582 continue                                                              295
10571 continue                                                              296
10572 continue                                                              296
      return                                                                297
      end                                                                   298
      subroutine lasso(rho,n,vv,s,thr,x,z,mm)                               299
      implicit double precision(a-h,o-z)                                    300
      double precision rho(n),vv(n,n),s(n),x(n),z(n)                        300
      integer mm(n)                                                         301
      call fatmul(2,n,vv,x,s,z,mm)                                          302
10590 continue                                                              302
10591 continue                                                              302
      dlx=0.0                                                               303
10600 do 10601 j=1,n                                                        303
      xj=x(j)                                                               303
      x(j)=0.0                                                              303
      t=s(j)+vv(j,j)*xj                                                     304
      if (abs(t)-rho(j).gt.0.0) x(j)=sign(abs(t)-rho(j),t)/vv(j,j)          305
      if(x(j).eq.xj)goto 10601                                              305
      del=x(j)-xj                                                           306
      dlx=max(dlx,abs(del))                                                 306
      s=s-del*vv(:,j)                                                       307
10601 continue                                                              308
10602 continue                                                              308
      if(dlx.lt.thr)goto 10592                                              308
      goto 10591                                                            309
10592 continue                                                              309
      return                                                                310
      end                                                                   311
      subroutine fatmul(it,n,vv,x,s,z,m)                                    312
      implicit double precision(a-h,o-z)                                    313
      parameter(fac=0.2)                                                    314
      double precision vv(n,n),x(n),s(n),z(n)                               314
      integer m(n)                                                          315
      l=0                                                                   316
10610 do 10611 j=1,n                                                        316
      if(x(j).eq.0.0)goto 10611                                             316
      l=l+1                                                                 316
      m(l)=j                                                                316
      z(l)=x(j)                                                             316
10611 continue                                                              317
10612 continue                                                              317
      if(l .le. int(fac*n))goto 10631                                       318
      if(it .ne. 1)goto 10651                                               318
      s=matmul(vv,x)                                                        318
      goto 10661                                                            318
10651 continue                                                              318
      s=s-matmul(x,vv)                                                      318
10661 continue                                                              319
10641 continue                                                              319
      goto 10621                                                            320
10631 if(it .ne. 1)goto 10671                                               320
10680 do 10681 j=1,n                                                        320
      s(j)=dot_product(vv(j,m(1:l)),z(1:l))                                 320
10681 continue                                                              320
10682 continue                                                              320
      goto 10691                                                            321
10671 continue                                                              321
10700 do 10701 j=1,n                                                        321
      s(j)=s(j)-dot_product(vv(m(1:l),j),z(1:l))                            321
10701 continue                                                              321
10702 continue                                                              321
10691 continue                                                              322
10621 continue                                                              322
      return                                                                323
      end                                                                   324
      subroutine inv(n,ww,xs,wwi)                                           325
      implicit double precision(a-h,o-z)                                    326
      double precision ww(n,n),xs(n-1,n),wwi(n,n)                           327
      nm1=n-1                                                               327
      xs=-xs                                                                328
      wwi(1,1)=1.0/(ww(1,1)+dot_product(xs(:,1),ww(2:n,1)))                 329
      wwi(2:n,1)=wwi(1,1)*xs(:,1)                                           330
      wwi(n,n)=1.0/(ww(n,n)+dot_product(xs(:,n),ww(1:nm1,n)))               331
      wwi(1:nm1,n)=wwi(n,n)*xs(:,n)                                         332
10710 do 10711 j=2,nm1                                                      332
      jm1=j-1                                                               332
      jp1=j+1                                                               333
      wwi(j,j)=1.0/(ww(j,j)+dot_product(xs(1:jm1,j),ww(1:jm1,j))  +dot_p    335 
     *roduct(xs(j:nm1,j),ww(jp1:n,j)))
      wwi(1:jm1,j)=wwi(j,j)*xs(1:jm1,j)                                     336
      wwi(jp1:n,j)=wwi(j,j)*xs(j:nm1,j)                                     337
10711 continue                                                              338
10712 continue                                                              338
      return                                                                339
      end                                                                   341
      subroutine glassopath(beta,what,jerrs,rholist,nrho,n,ss,rho,ia,itr    343 
     *, 		 ipen,thr,maxit,ww,wwi,niter,del,jerr)
      implicit double precision(a-h,o-z)                                    344
      integer nrho,n,jerrs(nrho)                                            345
      double precision rholist(nrho),beta(n,n,nrho),what(n,n,nrho)          346
      double precision ss(n,n),rho(n,n),ww(n,n),wwi(n,n)                    348
      is=0                                                                  350
10720 do 10721 j=1,n                                                        352
10730 do 10731 k=1,n                                                        353
      rho(j,k) = rholist(nrho)                                              354
10731 continue                                                              355
10732 continue                                                              355
10721 continue                                                              356
10722 continue                                                              356
      call glasso(n,ss,rho,ia,is,itr,ipen,thr,maxit,ww,wwi,niter,del,jer    357 
     *r)
      jerrs(1)=jerr                                                         359
10740 do 10741 j=1,n                                                        361
10750 do 10751 k=1,n                                                        362
      beta(j,k,nrho)=wwi(j,k)                                               363
      what(j,k,nrho)=ww(j,k)                                                364
10751 continue                                                              365
10752 continue                                                              365
10741 continue                                                              366
10742 continue                                                              366
      is=1                                                                  368
10760 do 10761 i=nrho,1,-1                                                  370
10770 do 10771 j=1,n                                                        372
10780 do 10781 k=1,n                                                        373
      rho(j,k)=rholist(i)                                                   374
10781 continue                                                              375
10782 continue                                                              375
10771 continue                                                              376
10772 continue                                                              376
      if(itr.gt.0) call dblepr('rho=', -1, rholist(i), 1)                   377
      itr2=itr                                                              378
      if(itr2.gt.0) itr2=itr-1                                              379
      call glasso(n,ss,rho,ia,is,itr2,ipen,thr,maxit,ww,wwi,niter,del,je    380 
     *rr)
      jerrs(i)=jerr                                                         382
10790 do 10791 j=1,n                                                        384
10800 do 10801 k=1,n                                                        385
      beta(j,k,i)=wwi(j,k)                                                  386
      what(j,k,i)=ww(j,k)                                                   387
10801 continue                                                              388
10802 continue                                                              388
10791 continue                                                              389
10792 continue                                                              389
10761 continue                                                              390
10762 continue                                                              390
      return                                                                391
      end                                                                   393
