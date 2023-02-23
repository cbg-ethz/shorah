#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <x86intrin.h>

int main(int argc, char * argv[])
{
    /* the following test is for emulated 32-bit on physical 64-bit */
    if (sizeof(unsigned long) != 8)
      abort ();
    volatile uint64_t a0 = UINT64_C(0x12389714dab0217);
    return _mm_popcnt_u64 (a0) == 25 ? EXIT_SUCCESS : EXIT_FAILURE;
}
