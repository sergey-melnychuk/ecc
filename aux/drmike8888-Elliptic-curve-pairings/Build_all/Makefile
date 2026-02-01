GMP_PREFIX=/opt/homebrew/opt/gmp
CFLAGS=-I$(GMP_PREFIX)/include -Ik12.header
LDFLAGS=-L$(GMP_PREFIX)/lib -L/usr/local/lib

OBJS = modulo.o poly.o eliptic.o poly_eliptic.o pairing.o signature.o base_protocols.o snarkbase.o

BINS = test_mod test_irrd pairing_sweep_alpha get_curve eliptic_test \
       weil_6_bit_pairing tate_6_bit_pairing quotient_group \
       signatures_11 signatures_11_keygen \
       base_curves_embed base_test base_curve_gen poly_exp_test \
       snark_qap snark_crs snark_proof snark_verify

all: $(BINS)
	@echo "All binaries built successfully."

test: $(BINS)
	@echo "=== Running tests ==="
	-@echo "\n--- test_mod ---" && ./test_mod
	-@echo "\n--- test_irrd ---" && ./test_irrd
	-@echo "\n--- poly_exp_test 5 ---" && ./poly_exp_test 5
	-@echo "\n--- eliptic_test 1847 23 42 ---" && ./eliptic_test 1847 23 42
	-@echo "\n--- weil_6_bit_pairing ---" && ./weil_6_bit_pairing
	-@echo "\n--- tate_6_bit_pairing ---" && ./tate_6_bit_pairing
	-@echo "\n--- quotient_group ---" && ./quotient_group
	-@echo "\n--- signatures_11_keygen ---" && ./signatures_11_keygen
	-@echo "\n--- signatures_11 ---" && ./signatures_11
	-@echo "\n--- base_test ---" && ./base_test
	-@echo "\n--- snark_qap ---" && ./snark_qap
	-@echo "\n--- snark_crs ---" && ./snark_crs
	-@echo "\n--- snark_proof ---" && ./snark_proof
	-@echo "\n--- snark_verify ---" && ./snark_verify
	@echo "\n=== All tests completed ==="

modulo.o:	modulo.c modulo.h
	gcc $(CFLAGS) -c -o modulo.o modulo.c

poly.o:  	poly.h poly.c
	gcc $(CFLAGS) -c -o poly.o poly.c

eliptic.o:	eliptic.c eliptic.h
	gcc $(CFLAGS) -c -o eliptic.o eliptic.c

poly_eliptic.o: poly_eliptic.h poly_eliptic.c
	gcc $(CFLAGS) -c -o poly_eliptic.o poly_eliptic.c

test_mod:	test_mod.c modulo.o poly.o
	gcc $(CFLAGS) -o test_mod test_mod.c modulo.o poly.o $(LDFLAGS) -lgmp

test_irrd:	test_irrd.c modulo.o poly.o eliptic.o poly_eliptic.o
	gcc $(CFLAGS) -o test_irrd test_irrd.c modulo.o poly.o eliptic.o poly_eliptic.o $(LDFLAGS) -lgmp

pairing_sweep_alpha:	pairing_sweep_alpha.c
	gcc $(CFLAGS) -o pairing_sweep_alpha pairing_sweep_alpha.c $(LDFLAGS) -lgmp -lm

get_curve:	get_curve.c modulo.o poly.o eliptic.o
	gcc $(CFLAGS) -o get_curve get_curve.c modulo.o poly.o eliptic.o $(LDFLAGS) -lgmp

eliptic_test:	eliptic_test.c modulo.o eliptic.o
	gcc $(CFLAGS) -o eliptic_test eliptic_test.c eliptic.o modulo.o $(LDFLAGS) -lgmp

pairing.o:	pairing.c pairing.h
	gcc $(CFLAGS) -c -o pairing.o pairing.c

weil_6_bit_pairing:	weil_6_bit_pairing.c modulo.o eliptic.o poly.o poly_eliptic.o pairing.o
	gcc $(CFLAGS) -o weil_6_bit_pairing weil_6_bit_pairing.c modulo.o poly.o eliptic.o poly_eliptic.o pairing.o $(LDFLAGS) -lgmp

tate_6_bit_pairing:	tate_6_bit_pairing.c modulo.o eliptic.o poly.o poly_eliptic.o pairing.o
	gcc $(CFLAGS) -o tate_6_bit_pairing tate_6_bit_pairing.c modulo.o poly.o eliptic.o poly_eliptic.o pairing.o $(LDFLAGS) -lgmp

quotient_group:	quotient_group.c modulo.o eliptic.o poly.o poly_eliptic.o pairing.o
	gcc $(CFLAGS) -o quotient_group quotient_group.c modulo.o poly.o eliptic.o poly_eliptic.o pairing.o $(LDFLAGS) -lgmp

signature.o: signature.c signature.h
	gcc $(CFLAGS) -c -o signature.o signature.c

signatures_11: signatures_11.c modulo.o eliptic.o poly.o poly_eliptic.o pairing.o signature.o
	gcc $(CFLAGS) -o signatures_11 signatures_11.c  modulo.o eliptic.o poly.o poly_eliptic.o pairing.o signature.o $(LDFLAGS) -lgmp -lk12

signatures_11_keygen: signatures_11_keygen.c modulo.o eliptic.o poly.o poly_eliptic.o pairing.o signature.o
	gcc $(CFLAGS) -o signatures_11_keygen signatures_11_keygen.c  modulo.o eliptic.o poly.o poly_eliptic.o pairing.o signature.o $(LDFLAGS) -lgmp -lk12

base_curves_embed:  base_curves_embed.c modulo.o
	gcc $(CFLAGS) -o base_curves_embed base_curves_embed.c modulo.o $(LDFLAGS) -lgmp

base_protocols.o: base_protocols.c base_protocols.h
	gcc $(CFLAGS) -c -o base_protocols.o base_protocols.c

base_test:  base_test.c base_protocols.o modulo.o eliptic.o
	gcc $(CFLAGS) -o base_test base_test.c base_protocols.o eliptic.o modulo.o $(LDFLAGS) -lgmp -lk12

base_curve_gen: base_curve_gen.c modulo.o eliptic.o
	gcc $(CFLAGS) -o base_curve_gen base_curve_gen.c eliptic.o modulo.o $(LDFLAGS) -lgmp

poly_exp_test:	poly_exp_test.c poly.o modulo.o
	gcc $(CFLAGS) -o poly_exp_test poly_exp_test.c poly.o modulo.o $(LDFLAGS) -lgmp

snarkbase.o: snarkbase.c modulo.h snarkbase.h
	gcc $(CFLAGS) -c -o snarkbase.o snarkbase.c

snark_qap: snark_qap.c snarkbase.o modulo.o eliptic.o poly.o poly_eliptic.o pairing.o signature.o
	gcc $(CFLAGS) -o snark_qap snark_qap.c snarkbase.o modulo.o eliptic.o poly.o poly_eliptic.o pairing.o signature.o $(LDFLAGS) -lgmp -lk12

snark_crs: snark_crs.c snarkbase.o modulo.o eliptic.o poly.o poly_eliptic.o pairing.o signature.o
	gcc $(CFLAGS) -o snark_crs snark_crs.c snarkbase.o modulo.o eliptic.o poly.o poly_eliptic.o pairing.o signature.o $(LDFLAGS) -lgmp -lk12

snark_proof: snark_proof.c snarkbase.o modulo.o eliptic.o poly.o poly_eliptic.o pairing.o signature.o
	gcc $(CFLAGS) -o snark_proof snark_proof.c snarkbase.o modulo.o eliptic.o poly.o poly_eliptic.o pairing.o signature.o $(LDFLAGS) -lgmp -lk12

snark_verify: snark_verify.c snarkbase.o modulo.o eliptic.o poly.o poly_eliptic.o pairing.o signature.o
	gcc $(CFLAGS) -o snark_verify snark_verify.c snarkbase.o modulo.o eliptic.o poly.o poly_eliptic.o pairing.o signature.o $(LDFLAGS) -lgmp -lk12

clean:
	rm -f $(OBJS) $(BINS)

.PHONY: all test clean
