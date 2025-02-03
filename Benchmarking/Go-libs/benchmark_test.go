package main

import (
	"testing"

	bls12_gnark "github.com/consensys/gnark-crypto/ecc/bls12-381"
	bls12_kilic "github.com/kilic/bls12-381"
)


func init() {
	// Load scalars for Benchmark
	for i := 0; i < 1000; i++ {
		scalars[i].SetString(BENCH_DATA[i], 10)
	}
}

// BenchmarkGnarkScalarMultiplication benchmarks the gnark library's scalar multiplication.
func BenchmarkGnarkScalarMultiplication(b *testing.B) {
	// Default Generator for BLS12-381
	var gj bls12_gnark.G1Jac
	gj.X.SetString("3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507")
	gj.Y.SetString("1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569")
	gj.Z.SetOne()

	// Reset the timer to exclude setup time
	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		// Use a scalar from the preloaded scalars array
		scalar := &scalars[i%1000]
		gj.ScalarMultiplication(&gj, scalar)
	}
}

// BenchmarkKilicScalarMultiplication benchmarks the Kilic library's scalar multiplication.
func BenchmarkKilicScalarMultiplication(b *testing.B) {
	// Default Generator for BLS12-381
	G1 := bls12_kilic.NewG1()
	g := bls12_kilic.G1One
	buf := G1.New()

	// Reset the timer to exclude setup time
	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		// Use a scalar from the preloaded scalars array
		scalar := &scalars[i%1000]
		G1.MulScalarBig(buf, &g, scalar)
	}
}