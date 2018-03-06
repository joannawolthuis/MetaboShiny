//The printspectrum tool prints out the spectrum (mz and intensity values) of a
//Thermo RAW File
//
//  Every line of the output is a peak registered by the mass spectrometer
//  characterized by an m/z value in Da and an intensity in the mass spectrometer's unit of abundance
package main

import (
	"../../unthermo"
	"flag"
	"fmt"
	"log"
)

func main() {
	var filename string

	//Parse arguments
	flag.StringVar(&filename, "raw", "small.RAW", "name of the RAW file")
	flag.Parse()

	//open RAW file
	file, err := unthermo.Open(filename)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Println(file.Scans)

	for i := uint64(1); i < file.Scans; i++ {
		scan := file.Scan(int(i))
		fmt.Println(scan.Polarity)
	}

	defer file.Close()
}
