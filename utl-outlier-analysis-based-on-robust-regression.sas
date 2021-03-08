Outlier analysis based on robust regression

Macro on end and in tools.

GitHub
https://tinyurl.com/nansw3z9
https://github.com/rogerjdeangelis/utl-outlier-analysis-based-on-robust-regression

VOODOO
https://github.com/rogerjdeangelis/voodoo

Tools
https://github.com/rogerjdeangelis/utl-macros-used-in-many-of-rogerjdeangelis-repositories

Inspired by
https://tinyurl.com/865xn2t6
https://communities.sas.com/t5/Statistical-Procedures/detect-outliers-using-macro-function/m-p/724302

The code is sensitive to the number of observatons and
provides an estimate of the expected number of outliers.
(random sample from normal and 1% tail)

I pulled this nmacro out of the voodoo package

*_                   _
(_)_ __  _ __  _   _| |_
| | '_ \| '_ \| | | | __|
| | | | | |_) | |_| | |_
|_|_| |_| .__/ \__,_|\__|
        |_|
;
data cars;
 set sashelp.cars(keep=msrp);
run;quit;

 Up to 40 obs WORK.CARS total obs=428

Obs     MSRP    HORSEPOWER    MPG_CITY    WEIGHT

  1    36945        265          17        4451
  2    23820        200          24        2778
  3    26990        200          22        3230
..

options ls=120 ps=32;
proc chart data=cars;
 hbar msrp;
run;quit;


    MSRP
  Midpoint                                                                         Freq
           |
 $15,000   |********************************************************                139
           |
 $30,000   |*********************************************************************   172
           |
 $45,000   |*****************************                                            72
           |
 $60,000   |*******                                                                  18
           |
 $75,000   |******                                                                   16
           |
 $90,000   |***                                                                       7
           |
$105,000   |                                                                          0
           |
$120,000   |*                                                                         2
           |
$135,000   |.                                                                         1
           |
$150,000   |                                                                          0
           |
$165,000   |                                                                          0
           |
$180,000   |                                                                          0
           |
$195,000   |.   -->looks like an outlier                                              1
           |
           ----+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+-
               10  20  30  40  50  60  70  80  90 100 110 120 130 140 150 160 170

*
 _ __  _ __ ___   ___ ___  ___ ___
| '_ \| '__/ _ \ / __/ _ \/ __/ __|
| |_) | | | (_) | (_|  __/\__ \__ \
| .__/|_|  \___/ \___\___||___/___/
|_|
;

%_vdo_outlyr(lib=work,mem=cars);

*            _               _
  ___  _   _| |_ _ __  _   _| |_
 / _ \| | | | __| '_ \| | | | __|
| (_) | |_| | |_| |_) | |_| | |_
 \___/ \__,_|\__| .__/ \__,_|\__|
                |_|
;

 
10 worst outliers up out of 27 outliers (expected_outliers=5.9212240204)
Robust Regression with 2.4615286299 * sigma cuttoff

Obs    OBS    SIGMAS       MSRP

  1    335    10.9085    $192,465

  2    263     6.6545    $128,420
  3    272     6.5383    $126,670
  4    271     6.2128    $121,770
  5    262     4.4227     $94,820
  6    270     4.1371     $90,520
  7      7     4.0870     $89,765
  8    201     3.9030     $86,995
  9    269     3.9013     $86,970
 10     21     3.7439     $84,600

*
 _ __ ___   __ _  ___ _ __ ___
| '_ ` _ \ / _` |/ __| '__/ _ \
| | | | | | (_| | (__| | | (_) |
|_| |_| |_|\__,_|\___|_|  \___/

;

%macro _vdo_outlyr(lib=&libname,mem=&data);

  /*
    %let libname =sashelp;
    %let data=bweight;
    %let var=weight;
  */

 %LOCAL
         UIN
         RC
         UI
         UDSID
         UOBS
         UVARS
         ULBL
         VARNAM
         VARTYPE
         VARLABEL
         VARLEN
         OBS
         EXPECTED OUTLIERS
         OBS
         NOBS
         DOTHEREST
;

 %let uin=&lib..&mem;

 %put %sysfunc(ifc(%sysevalf(%superq(uin)=,boolean),**** Please Provide SAS dataset ****,));

 %if %sysfunc(ifc(%sysevalf(%superq(uin)=,boolean),1,0)) eq 0 %then %do;

    %local nobs cutoff expected_outliers dotherest;

    * grubs test is better but it is too computationally expensive
      It recomputes thowing ot the worst outlier;

    *  I do have a couple of grub test macros in R and SAS;

    * arguments have already been checked by the call driver;

    * estimate the best cuttoff in terms of the sigma on standardized values;

    * number of obs in source data;
    proc sql noprint; select count(*) into :nobs from &uin;quit;

    %put &=nobs;

    * get the best cutoff in terms of N * sigma and the expected outliers;
    data _null;
      * (e^x)^2) -> sqrt(log(x));
      cutoff=sqrt(log(&nobs));
      call symputx('cutoff',cutoff);
      put cutoff=;
      expected_outliers=2 * &nobs *(1-cdf('normal',cutoff));
      call symputx('expected_outliers',expected_outliers);
      put expected_outliers=;
    run;quit;

    %LET UDSID = %SYSFUNC( OPEN ( &uin., I ) );

   /*-------------------------------------*\
   ! GET THE NUMBER OF COLUMNS FOR LOOP    !
   \*-------------------------------------*/

   %LET UOBS = %SYSFUNC(ATTRN(&UDSID,NLOBS));

   %LET UVARS = %SYSFUNC(ATTRN(&UDSID,NVARS));

   %LET ULBL = %SYSFUNC(ATTRC(&UDSID,LABEL));

   %LET UMAX = 0;

   %DO UI = 1 %TO &UVARS;


     /*-------------------------------------*\
     !  GET ATTRIBUTES                       !
     \*-------------------------------------*/

     /* %let ui=1;  */

     %LET UVARNAM = %SYSFUNC ( VARNAME  ( &UDSID, &UI  ) );
     %LET UVARTYP = %SYSFUNC ( VARTYPE  ( &UDSID, &UI  ) );
     %LET UVARLBL = %SYSFUNC ( VARLABEL ( &UDSID, &UI  ) );
     %LET UVARLEN = %SYSFUNC ( VARLEN   ( &UDSID, &UI  ) );

     %IF &UVARTYP EQ %QUPCASE(N) %THEN %DO;

       %let var=&UVARNAM;

       * robust regression - only interested in outliers - leverage values might also be interesting;
       ods exclude all;
       ods output diagnostics=__vvdag;
       proc robustreg data=&uin method=MM;
       model &var = /diagnostics  cutoff=&cutoff /* stadardized sigma &cutoff * sigma */;
       run;
       ods select all;
       ods listing;

       * number of potential outliers;
       proc sql noprint;select count(*) into :obs from __vvdag;quit;

       * do we have more than the expected outliers;
       %let dotherest=%sysfunc(ifn(&obs > &expected_outliers,1,0));

       %put &=dotherest;

           * set up for sort by abs value;
           data __vvdagabs;
            set __vvdag;
            rresidual=abs(rresidual);
           run;quit;

           /* sort by abs value so we can remove the lower expected_outliers */
           proc sort data=__vvdagabs out=__vvdagsrt noequals;
           by rresidual;
           run;quit;

           /* drop the lower values */
           data __vvdagsel;
            set __vvdagsrt(firstobs=%sysfunc(int(&expected_outliers))); * remove expected outliers;
           run;quit;

           * go back to full data and get orginal value;
           * get bad values;
           data __vvdagget;
             do until (dne);
               set __vvdagsel(drop=outlier) end=dne;
               rec=obs;
               set &uin.(keep=&var)  point=rec;
               output;
             end;
             stop;
           run;quit;

           proc sort data=__vvdagget  out=__vvdagfin noequals;
           by descending rresidual;
           run;quit;

           title1 ' ';title2 ' ';title3 ' ' ;
           TITLE4 "10 worst outliers up out of &obs outliers (expected_outliers=&expected_outliers)";
           TITLE5 "Robust Regression with &cutoff * sigma cuttoff and removal of expected outliers?";

           proc print data=__vvdagfin(obs=10 rename=rresidual=sigmas) width=min;
           run;quit;
           title;
   proc datasets lib=work nolist;
   delete __vvdag:;
   run;quit;
   %end;
  %end;
%end;

%mend _vdo_outlyr;

/*

data car;
 set sashelp.cars(keep=msrp horsepower mpg_city weight);
run;quit;

%_vdo_outlyr(lib=work,mem=cars);

*/
