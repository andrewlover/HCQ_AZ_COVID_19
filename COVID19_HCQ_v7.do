pwd

*********************************************************************
* Project:  	COVID-19 HCQ				
* This file:	Survival analysis
* When/Who:	 	AAL 20 Mar 2020
*********************************************************************

cd "~/Documents"

version 16
clear all
macro drop _all    
rcall setpath "/usr/local/bin/R"
* set scheme plottig

***** Data source: Preprint
* Gautret et al (2020) IJAA, DOI: 10.1016/j.ijantimicag.2020.105949

******************************************************************
* Section I:  LTF additions
******************************************************************

/*
** from text:
Six hydroxychloroquine-treated patients were LTF
three patients were transferred to intensive care unit, including:
 
one transferred on day2 post-inclusion who was PCR-positive on day1, 
one transferred on day3 post-inclusion who was PCR-positive on days1-2 and 
one transferred on day4 post-inclusion who was PCR- positive on day1 and day3; 

one patient died on day3 post inclusion and was PCR-negative on day2; 
one patient decided to leave the hospital on day3 post-inclusion and was PCR-negative on days1-2
one patient stopped the treatment on day3 post-inclusion because of nausea and was PCR-positive on days 1-2-3

37. to ICU on day 2, PCR-pos day 1
38. to ICU on day 3, PCR-pos day 1-2
39. to ICU on day 4, PCR pos day 1/3
40. death on day 3, PCR-neg day 2
41. LTF on day 3, pcr neg day 2
42. withd day 3, pcr + 1-3

*/

input patient  str3 cq_treat str8 d0 str8 d1 str8 d2 str8 d3 str8 d4 d_censor str20 note
*		HCQ		D0		D1		D2		D3		D4

37 		Yes 	"POS"	"POS" 	"ND" 	"ND" 	"ND" 2 	"ICU"
38 		Yes 	"POS" 	"POS" 	"POS" 	"ND" 	"ND" 3 	"ICU"
39 		Yes 	"POS" 	"POS" 	"ND" 	"POS" 	"ND" 4 	"ICU"

40 		Yes 	"POS" 	"POS" 	"NEG" 	"ND" 	"ND" 3 	"death"
41 		Yes 	"POS" 	"POS" 	"NEG" 	"ND" 	"ND" 3 	"LTF"
42 		Yes 	"POS" 	"POS" 	"POS" 	"POS" 	"ND" 3 	"withdrawn"
end

save "~/Documents/_HCQ_LTF_data_v1.dta", replace

******************************************************************
* Section II:  data cleaning
******************************************************************

import delimited "~/Desktop/COVID_19_HCQ/1_DATA/tabula-HCQ__v2y.csv", clear // extracted via Tabula

append using "~/Documents/_HCQ_LTF_data_v1.dta"

drop v5 v9 
replace sex = "unknown",, if sex == ""
sencode sex, replace

* recode to be 1 == virus-free 0 == detectable viremia
**  negative PCR (CT value >= 35) -----> 1

foreach var of varlist d0 - d6 {
	list patient `var'
}	
*

// save "recode_test.dta", replace
 
foreach var of varlist d0- d6 {
	list d0 - d6
	replace `var' = "999",, if `var' == "NEG" // standin for '1'
	list d0 - d6
	replace `var' = "0",, if `var' == "POS"
	replace `var' = ".",, if `var' == "ND"
	replace `var' = ".",, if `var' == "NF"
	list d0 - d6
	destring `var', replace
	replace `var' = 1,, if (`var' >= 35) & ((`var' != 999) & (`var' != .))
	replace `var' = 0,, if (`var' < 35) & ((`var' != 999) & (`var' != .))
	replace `var' = 1,, if (`var' == 999)
	list d0 - d6, sep(5)
	* sleep 5000
}
*

tab d5, miss

gen arm = ""

replace arm = "control" if cq == "No" & azith == "No"
replace arm = "HCQ" if cq == "Yes" 
tab arm

gen arm2 = ""

replace arm2 = "control" if cq == "No" & azith == "No"
replace arm2 = "HCQ" if cq == "Yes" 
replace arm2 = "HCQ+AZ" if cq == "Yes" & azith == "Yes"

tab arm2

gsort -arm
sencode arm, replace

gsort arm2
sencode arm2, replace

split hcq_serum, parse(" ") gen (hcq_s)
replace hcq_s1 = "",, if hcq_s1 == "-"

destring hcq_s1, replace

charlist hcq_s2
return list
display r(ascii)
replace hcq_s2 = subinstr(hcq_s2, "`=char(40)'", "",.) 
replace hcq_s2 = subinstr(hcq_s2, "`=char(41)'", "",.) 
replace hcq_s2 = subinstr(hcq_s2, "`=char(68)'", "",.) 

listsome hcq_s2, random
destring hcq_s2, replace

tab clinicalstatus

replace clinicalstatus = "Asymptomatic",, if clinicalstatus == "Asvmotomatic" | clinicalstatus == "Asvmntomatic"
replace clinicalstatus = "LTRI",, if clinicalstatus == "LRTI"

sencode clinicalstatus, replace

******************************************************************
* Section II:  Primary Outcome
******************************************************************

* Gautret et al:
* "The primary endpoint was virological clearance at day-6 post-inclusion."

tab d6, miss
binreg d6 i.arm i.sex age, rr
outreg2 using "HCQ_2020_T1.tex", ///
	replace sideway stats(coef ci pval) noaster dec(3) label tex eform
	

// binreg d6 i.arm2 i.sex age, or // perfect prdiction

logistic d6 i.arm
estat ic
tjur2

logistic d6 i.arm ib2.clinical // not supported by BIC
estat ic
tjur2

// binreg d6 i.arm ib2.clinical i.sex age, rr ml diff

firthlogit d6 i.arm i.sex age, or
firthfit

firthlogit d6 i.arm ib1.clinical i.sex age, or // not supported
firthfit

firthlogit d6 ib3.arm2 i.sex age, or
firthfit

outreg2 using "HCQ_2020_T2.tex", ///
	replace sideway stats(coef ci pval) noaster dec(3) label tex eform
	
	tab d6,  miss
		

******************************************************************
* Section III:  Secondary outcome
******************************************************************

* Gautret et al:
* "Secondary endpoint was virological clearanc e primary endpoint was virological clearance at day-6 post-inclusion."

tab d0 d1, miss
tab d0 d6, miss

*****
reshape long d, i(patient) j(day)
*****

tab d, miss

rename d vfree

sort patient day, stable

tsset patient day

// ssc install tsspell
tsspell, cond(vfree > 0 & vfree < .)

listsome p day vfree _*

gen failtime1 = . // first PCR(-) day
replace failtime1 = day,, if _seq == 1 & _spell == 1

foreach var of varlist _seq - _end {
	rename `var' `var'1

}
*

gen seq2 = . 
by patient: replace seq2 = 1,, if (vfree == 1) & (vfree[_n-1] == 1) & (vfree[_n-2] == 0)

listsome p day vfree seq2,, if patient == 1

gen failtime2 = . // two sequential PCR(-) days
replace failtime2 = day,, if seq2 == 1

missings dropvars, force

order patient vfree* fail* age arm

* need to drop extra failtimes

bysort patient: egen t1 = min(failtime1)
bysort patient: egen t2 = min(failtime2)

rename patient id

by id: tab vfree

drop fail* _seq seq _spell _end

******************************************************************
* Section IV:  st etc
******************************************************************
*** needed to get sqindexplot axis ticks
gen day2 = day + 1
tab day2

tostring day2, replace

tab day2

// tostring day2, replace
// replace day2 = "D0",, if day == 0
// tab day2

label def day2 1 "D0" 2 "1" 3 "2" 4 "3" 5 "4" 6 "5" 7 "6"
sencode day2, label(day2) replace

tab day2, nolab

tab day day2

tab day2, nolab
sqset vfree id day

* 1-16 control; 17-36 HCQ; 37-42 LTF

tab id

sort id, stable

tostring id, gen(id2)
sencode id2, replace

tostring vfree, replace
replace vfree = "missing" if vfree == "."

sencode vfree, replace

sqindexplot, color(red blue green black) gapinclude order(id)

* 1-16 control; 17-36 HCQ; 37-42 LTF

sqindexplot, scheme(plottig) gapinclude order(id) color("sea" "plg3" gs10) rbar ///
		xtitle("Study day", margin(vsmall) size(4.5)) ///
		xlabel(0 "D0" 1 "1" 2 "2" 3 "3" 4 "4" 5 "5" 6 "6", valuelabel labsize(4)) ///
		ytitle("Patient ID", size(4.5)) ///
		ylabel(1(3)42, angle(h) labsize(4)) ///
		text(4 -1.5 "{bf:Control}", color(maroon) place(se) just(left)) ///
		text(30 -1.5 "{bf:HCQ}", color(maroon) place(se) just(left)) ///
		text(40 -1.5 "{bf:LTF}", color(maroon) place(se) just(left)) ///
		legend(label(1 "PCR(+) Ct < 35") label(3 "PCR(-) Ct >=35") label(5 "Missing data/LTF")) ///
		legend(order(1 3 5) col(3) size(4) position(11)) ///
		title("") ///
		yline(16 36, lwidth(1.2) lcolor(maroon) lpattern(solid)) ///
		note("Reanalysis of data from Gautret {it:et al.} (2020) IJAA, DOI: 10.1016/j.ijantimicag.2020.105949" ///
			/*"@AndrewALover"*/, size(3) pos(12))

graph export "Gautret_HCQ_2020_Seq_v1.pdf", replace

drop day2

reshape wide vfree, i(id) j(day)

list id v* d_censor t1 t2

gen t_fail = 6
replace t_fail = t1,, if t1 != .
replace t_fail = d_censor,, if d_censor !=.
tab t_fail

gen fail1 = 0
replace fail1 = 1,, if t1 !=.
tab fail1

stset t_fail, fail(fail1)

list t1 t_fail fail1 d_censor _*

sts test arm2, logrank

sts graph, by(arm2) ///
		xtitle("Day", margin(tiny) size(large)) ///
		xlabel(0 "D0" 1 "1" 2 "2" 3 "3" 4 "4" 5 "5" 6 "6", valuelabel labsize(4)) ///
		legend(label (3 "Control") label (1 "HCQ") label (2 "HCQ-AZ")) ///
		legend(order(1 2 3)) scheme (plottig) /// 
		xtitle("Days", margin(vsmall) size(5)) ///
		ytitle("Proportion with undetectable" "viremia via PCR", size(5)) ///
		xlabel(, angle(h) labsize(4)) ///
		ylabel(0(0.2)1.0, angle(h) format(%3.1f) labsize(5)) ///
		legend(order(3 2 1) row(3) ring(0) size(4) position(8)) ///
		title("") ///
		plot1opts(lpattern(solid) lwidth(2) lcolor(navy)) ///
		plot2opts(lpattern(solid) lwidth(2) lcolor(gs8)) ///
		plot3opts(lpattern(solid) lwidth(2) lcolor("166 86 40")) ///
		note("Reanalysis of data from Gautret {it:et al.} (2020) IJAA, DOI:10.1016/j.ijantimicag.2020.105949" ///
			/*"@AndrewALover"*/, size(3.2) pos(12))

graph export "Gautret_HCQ_2020_KM_v1.pdf", replace

stcox i.arm age i.sex ib2.clinical
stcox ib3.arm2 age i.sex ib2.clinical

stcox ib3.arm2 // age i.sex ib2.clinical
 
stcox ib3.arm2 age i.sex // ib2.clinical

outreg2 using "HCQ_2020_T3.tex", ///
	replace sideway stats(coef ci pval) noaster dec(3) label tex eform

exit

stcurve, surv at1(arm2=1) at2(arm2=2) at3(arm2=3) ci

stcurve, surv at1(arm2=1) at2(arm2=2) at3(arm2=3) lwidth(1.5 1.5 1.5) ///
		legend(label (1 "Control") label (2 "HCQ") label (3 "HCQ-AZ")) ///
		legend(order(1 2 3)) scheme (plottig) /// 
		xtitle("Days", margin(vsmall) size(5)) ///
		ytitle("Proportion with undetectable" "viremia via PCR", size(5)) ///
		xlabel(, angle(h) labsize(4)) ///
		ylabel(0(0.2)1.0, angle(h) format(%3.1f) labsize(5)) ///
		legend(order(1 2 3) col(1) stack ring(0) size(4) position(8)) ///
		title("") ///
		note("Data from Gautret {it:et al.} (2020) IJAA, DOI : 10.1016/j.ijantimicag.2020.105949", size(3) pos(12))

graph export "Gautret_HCQ_2020_F1_v1.pdf", replace

exit


******************************************************************
* Section III:  
******************************************************************



******************************************************************
* Section IV:  
******************************************************************



******************************************************************
* Section V:  
******************************************************************



******************************************************************
* Section VII:  
******************************************************************



******************************************************************
* Section VIII:  
******************************************************************


exit

