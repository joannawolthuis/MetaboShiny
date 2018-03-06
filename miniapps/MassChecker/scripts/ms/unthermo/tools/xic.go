/*The XIC tool prints mass chromatograms for a specific m/z.

  For the m/z given, it prints the peak with highest intensity in interval
  [mz-tol ppm,mz+tol ppm] for every MS-1 scan.

  Every line contains the retention time and intensity of a found peak

  Example:
      xic -mz 361.1466 -tol 2.5 -raw rawfile.raw

  Output:
      0.003496666666666667 10500.583
      0.015028333333333333 11793.04
      0.03391333333333333 10178.598
      0.05393333333333334 10671.821
      0.07350833333333334 11572.251
*/
package main

import (
	"bitbucket.org/proteinspector/ms"
	"bitbucket.org/proteinspector/ms/unthermo"
	"flag"
	"fmt"
	"log"
	"sort"
)

//mz is the m/z for the XIC
var mz float64

//tol is the tolerance in ppm
var tol float64

//fileName is the file name of the raw file
var fileName string

func init() {
	flag.Float64Var(&mz, "mz", 810.41547, "m/z to filter on, this flag may be specified multiple times")
	flag.Float64Var(&tol, "tol", 2.5, "allowed m/z tolerance in ppm, can be used with -mz")
	flag.StringVar(&fileName, "raw", "small.RAW", "name of the subject RAW file")
	flag.Parse()
}

/*
  XIC peaks get extracted out of each MS1 Scan read by the unthermo library
*/
func main() {
	file, err := unthermo.Open(fileName)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	for i := 1; i <= file.NScans(); i++ {
		printXICpeak(file.Scan(i), mz, tol)
	}
}

//printXICpeaks outputs mz, scan time and intensity of the highest MS1 peak
//with mz within tolerance around the supplied mz.
func printXICpeak(scan ms.Scan, mz float64, tol float64) {
	if scan.MSLevel == 1 {
		filteredSpectrum := mzFilter(scan.Spectrum(), mz, tol)
		if len(filteredSpectrum) > 0 {
			fmt.Println(scan.Time, maxPeak(filteredSpectrum).I)
		}
	}
}

//mzFilter outputs the spectrum within tol ppm around the supplied mz
func mzFilter(spectrum ms.Spectrum, mz float64, tol float64) ms.Spectrum {
	return mzIntervalFilter(spectrum, mz-10e-6*tol*mz, mz+10e-6*tol*mz)
}

//mzIntervalFilter filters the spectrum for mz's within the interval [min,max)
//including minMz and excluding maxMZ
func mzIntervalFilter(spectrum ms.Spectrum, minMz float64, maxMz float64) ms.Spectrum {
	//A spectrum is sorted by m/z so we can do binary search for two
	//endpoint peaks and get the peaks between them.
	lowi := sort.Search(len(spectrum), func(i int) bool { return spectrum[i].Mz >= minMz })
	highi := sort.Search(len(spectrum)-lowi, func(i int) bool { return spectrum[i+lowi].Mz >= maxMz })

	return spectrum[lowi : highi+lowi]
}

//maxPeak returns the maximally intense peak within the supplied spectrum
func maxPeak(spectrum ms.Spectrum) (maxIn ms.Peak) {
	for _, peak := range spectrum {
		if peak.I >= maxIn.I {
			maxIn = peak
		}
	}
	return
}
