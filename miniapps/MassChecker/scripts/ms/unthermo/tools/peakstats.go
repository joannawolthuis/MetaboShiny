/*The peakstats tool outputs a few data about the peaks of supplied ions:
  - Mass, Time at Maximum, Maximal intensity of ions found in the LC/MS map
  - Full width at half maximum of this maximal peak
*/
package main

import (
	"bitbucket.org/proteinspector/ms"
	"bitbucket.org/proteinspector/ms/unthermo"
	"flag"
	"fmt"
	"github.com/pkelchte/spline"
	"log"
	"sort"
)

//tol is the tolerance (in ppm) for m/z peaks
const tol = 2.5

type TimedPeak struct {
	ms.Peak
	Time float64
}

func main() {
	//ions are the m/z that will be searched in the spectra
	//var ions = []float64{495.78700, 424.25560, 507.81340, 461.74760, 740.40170, 820.47250, 682.34770} //BSA
	var ions = []float64{363.67450, 362.22910, 367.21590, 550.76660, 643.85824, 878.47842, 789.90439} //Enolase

	var fileName string
	flag.StringVar(&fileName, "raw", "small.RAW", "name of the subject RAW file")
	flag.Parse()

	file, err := unthermo.Open(fileName)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	xicmap := xics(file, ions)
	timeresolution := guessMsOneInterval(file) / 10

	//loop over found ions in order
	var keys []float64
	for k := range xicmap {
		keys = append(keys, k)
	}
	sort.Float64s(keys)

	for _, mz := range keys {
		maxPeak, fwhm := maxChromFeature(xicmap[mz], timeresolution)
		fmt.Println(mz, maxPeak.Time, maxPeak.I, fwhm)
	}
}

//xics returns a map of slices of extracted ion chromatograms
func xics(file unthermo.File, ions []float64) map[float64][]TimedPeak {
	xicmap := make(map[float64][]TimedPeak, len(ions))

	for i := 1; i <= file.NScans(); i++ {
		scan := file.Scan(i)

		if scan.MSLevel == 1 {
			spectrum := scan.Spectrum()
			for _, ion := range ions {
				filteredSpectrum := mzFilter(spectrum, ion, tol)
				if len(filteredSpectrum) > 0 {
					xicmap[ion] = append(xicmap[ion], TimedPeak{maxPeak(filteredSpectrum), scan.Time})
				}
			}
		}
	}
	return xicmap
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

//guessMsOneInterval returns a guess for the interval between MS1 scans.
//Usually, aqcuisitions start with a few MS1 scans after each other,
//the minimum time between MS1 scans is then the time between the first two
func guessMsOneInterval(file unthermo.File) float64 {
	var timeOne float64 = 0
	var i int = 1
	for ; timeOne == 0; i++ {
		scan := file.Scan(i)
		if scan.MSLevel == 1 {
			timeOne = scan.Time
		}
	}
	timeTwo := timeOne
	for ; timeTwo == timeOne; i++ {
		scan := file.Scan(i)
		if scan.MSLevel == 1 {
			timeTwo = scan.Time
		}
	}
	return timeTwo - timeOne
}

//maxChromFeature outputs the maximal peak in an ion chromatogram, along
//with the full-width half max.
//The resolution is the time step size when searching for the half max
func maxChromFeature(xic []TimedPeak, resolution float64) (max TimedPeak, fwhm float64) {
	//look for the maximum peak in the xic
	for _, peak := range xic {
		if peak.I >= max.I {
			max = peak
		}
	}

	//create spline for having a mathematical function when searching for the fwhm.
	//it covers the whole xic for code conciseness, there is not much speed
	//gain when only a region around max is interpolated.
	s := spline.Spline{}
	X := make([]float64, len(xic))
	Y := make([]float64, len(xic))
	for i := range xic {
		X[i] = xic[i].Time
		Y[i] = float64(xic[i].I)
	}
	s.Set_points(X, Y, true)

	//find the time points when the xic goes below half the maximum
	right := max.Time
	for ; s.Operate(right) > float64(max.I/2); right += resolution {
	}
	left := max.Time
	for ; s.Operate(left) > float64(max.I/2); left -= resolution {
	}

	return max, right - left
}
