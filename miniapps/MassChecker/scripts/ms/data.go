//Package ms is a library for mass spectrometry data
package ms

//Peak represents an ion peak
type Peak struct {
	Mz float64
	I  float32
}

//A Spectrum is a collection of peaks
type Spectrum []Peak

//Scan represents the peak acquisition event of the mass spectrometer
type Scan struct {
	Analyzer Analyzer
	MSLevel  uint8
	//Spectrum is a function forcing the read of a spectrum,
	//which is "delayed" for efficiency reasons. If it was not delayed
	//and Spectrum were a data structure, it would always have to
	//be read, which is very expensive. Now if only another property of
	//Scan (cheaper to obtain) is requested, resources are saved.
	Spectrum func(centroided ...bool) Spectrum
	//PrecursorMzs is only filled with mz values at MSx scans.
	PrecursorMzs []float64
	Time         float64
	Polarity	string
	Preamble	[132]uint8
}

//Analyzer is the mass analyzer
type Analyzer int

//The analyzer types are documented in literature
const (
	ITMS Analyzer = iota
	TQMS
	SQMS
	TOFMS
	FTMS
	Sector
	Undefined
)

//Spectrum implements sort.Interface for []Peak based on m/z

func (a Spectrum) Len() int           { return len(a) }
func (a Spectrum) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
func (a Spectrum) Less(i, j int) bool { return a[i].Mz < a[j].Mz }
