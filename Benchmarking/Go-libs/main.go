//+build !amd64 generic
package main

import (
	"fmt"
	"math/big"
	"time"

	bls12_gnark "github.com/consensys/gnark-crypto/ecc/bls12-381"
	bls12_kilic "github.com/kilic/bls12-381"
	"github.com/dterei/gotsc"
)

var scalars [1000]big.Int

func main() {
	// Load scalars for Benchmark
	for i := 0; i < 1000; i++ {
		scalars[i].SetString(BENCH_DATA[i], 10)
	}

	// Benchmarking gnark-library scalar multiplication
	// Default Generator for BLS12-381
	var gj bls12_gnark.G1Jac
	gj.X.SetString("3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507")
	gj.Y.SetString("1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569")
	gj.Z.SetOne()
	var totalDuration time.Duration
	totalDuration = 0
	avgc := 0
	for i := 0; i < 1000; i++ {
		start := time.Now()
		startc := gotsc.BenchStart()
		gj.ScalarMultiplication(&gj, &scalars[i])		
		totalDuration += time.Since(start)
		endc := gotsc.BenchEnd()
		avgc += int(endc - startc)
	}
	averageDuration := totalDuration / 1000
	averageCycle := avgc / 1000
	fmt.Println("Average Duration:", averageDuration)
	fmt.Println("Average Cycle:", averageCycle)
	
	
	// Benchmarking Kilic-library scalar multiplication
	// Default Generator for BLS12-381
	G1 := bls12_kilic.NewG1()
	g := bls12_kilic.G1One
	b := G1.New()
	totalDuration = 0
	avgc = 0
	for i := 0; i < 1000; i++ {
		start := time.Now()
		startc := gotsc.BenchStart()
		G1.MulScalarBig(b, &g, &scalars[i])
		totalDuration += time.Since(start)
		endc := gotsc.BenchEnd()
		avgc += int(endc - startc)
	}
	averageDuration = totalDuration / 1000
	averageCycle = avgc / 1000
	fmt.Println("Average Duration:", averageDuration)
	fmt.Println("Average Cycle:", averageCycle)

}
