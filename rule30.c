#include <stdio.h>

//Code for experimenting with rule 30

unsigned int rule30_basic()
{
  static unsigned int x = 0x00010000; //32-bit seed
  x = (x>>1) ^ (x | x<<1);
  return x;
}

//Rule 30 generator, wraps lsb and msb around.
unsigned int rule30_wrap()
{
  static unsigned int x = 0x00010000; //32-bit seed
  x = ((x>>1)^(x<<31)) ^ (x | ((x<<1)^(x>>31)));
  return x;
}

void printer32(unsigned int x)
{
  for(int j=32; j--;)
  {
    printf("%c", x>>j & 1 ? '1' : '-');
  }
  printf("\n");
}

void main()
{
  unsigned int x;

  for(int i=0; i<100; i++)
  {
    x = rule30_wrap();
    printf("%u\n", x);
  }
  for(int i=0; i<100; i++)
  {
    x = rule30_basic();
    printer32(x);
  }
}
