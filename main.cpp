/*
 *  Shannon-Fano coding algorithm
 *  by Sergey Tikhonov (st@haqu.net)
 * 
 *  Usage: shannon [OPTIONS] input [output]
 *    The default action is to encode input file.
 *    -d  Decode file.
 * 
 *  Examples:
 *    shannon input.txt
 *    shannon input.txt encoded.txt
 *    shannon -d encoded.txt
 */

#ifdef _WIN32
#define _CRT_SECURE_NO_DEPRECATE
#define NL "\r\n"
#else
#define NL "\n"
#endif

#include <string>
#include <vector>
#include <map>
#include <cassert>

using namespace std;

struct pnode
{
  char ch; // char
  float p; // probability
};

static int pnode_compare( const void *elem1, const void *elem2 )
{
  const pnode a = *(pnode*)elem1;
  const pnode b = *(pnode*)elem2;
  if( a.p < b.p ) return 1; // 1 - less (reverse for decreasing sort)
  else if( a.p > b.p ) return -1;
  return 0;
}

class Coder
{
private:
  int tsize; // table size (number of chars)
  pnode *ptable; // table of probabilities
  map<char, string> codes; // codeword for each char

public:
  void Encode( const char *inputFilename, const char *outputFilename )
  {
    map<char, int> freqs; // frequency for each char from input text
    int i;

    //  Opening input file
    //
    FILE *inputFile;
    inputFile = fopen( inputFilename, "r" );
    assert( inputFile );

    //  Counting chars
    //
    char ch; // char
    unsigned total = 0;
    while( fscanf( inputFile, "%c", &ch ) != EOF )
    {
      freqs[ch]++;
      total++;
    }
    tsize = (int)freqs.size();

    //  Building decreasing freqs table
    //
    ptable = new pnode[tsize];
    assert( ptable );
    float ftot = float(total);
    map<char, int>::iterator fi;
    for( fi=freqs.begin(),i=0; fi!=freqs.end(); ++fi,++i )
    {
      ptable[i].ch = (*fi).first;
      ptable[i].p = float((*fi).second) / ftot;
    }
    qsort( ptable, tsize, sizeof(pnode), pnode_compare );

    //  Encoding
    //
    EncShannon( 0, tsize-1 );

    //  Opening output file
    //
    FILE *outputFile;
    outputFile = fopen( outputFilename, "wb" );
    assert( outputFile );

    //  Outputing ptable and codes
    //
    printf( "%i"NL, tsize );
    fprintf( outputFile, "%i"NL, tsize );
    for( i=0; i<tsize; i++ )
    {
      printf("%c\t%f\t%s"NL, ptable[i].ch, ptable[i].p, codes[ptable[i].ch].c_str() );
      fprintf(outputFile, "%c\t%f\t%s"NL, ptable[i].ch, ptable[i].p, codes[ptable[i].ch].c_str() );
    }

    //  Outputing encoded text
    //
    fseek( inputFile, SEEK_SET, 0 );
    printf(NL);
    fprintf(outputFile, NL);
    while( fscanf( inputFile, "%c", &ch ) != EOF )
    {
      printf("%s", codes[ch].c_str());
      fprintf(outputFile, "%s", codes[ch].c_str());
    }
    printf(NL);

    //  Cleaning
    //
    codes.clear();
    delete[] ptable;

    //  Closing files
    //
    fclose( outputFile );
    fclose( inputFile );
  }
  
  void Decode( const char *inputFilename, const char *outputFilename )
  {
    //  Opening input file
    //
    FILE *inputFile;
    inputFile = fopen( inputFilename, "r" );
    assert( inputFile );

    //  Loading codes
    //
    fscanf( inputFile, "%i", &tsize );
    char ch, code[128];
    float p;
    int i;
    fgetc( inputFile ); // skip end line
    for( i=0; i<tsize; i++ )
    {
      ch = fgetc(inputFile);
      fscanf( inputFile, "%f %s", &p, code );
      codes[ch] = code;
      fgetc(inputFile); // skip end line
    }
    fgetc(inputFile); // skip end line

    //  Opening output file
    //
    FILE *outputFile;
    outputFile = fopen( outputFilename, "w" );
    assert( outputFile );

    //  Decoding and outputing to file
    //
    string accum = "";
    map<char, string>::iterator ci;
    while((ch = fgetc(inputFile)) != EOF)
    {
      accum += ch;
      for( ci=codes.begin(); ci!=codes.end(); ++ci )
        if( !strcmp( (*ci).second.c_str(), accum.c_str() ) )
        {
          accum = "";
          printf( "%c", (*ci).first );
          fprintf( outputFile, "%c", (*ci).first );
        }
    }
    printf(NL);

    //  Cleaning
    //
    fclose( outputFile );
    fclose( inputFile );
  }

private:
  void EncShannon( int li, int ri )
  {
    int i, isp;
    float p;
    float pfull;
    float phalf;

    if( li == ri )
    {
      return;
    }
    else if( ri - li == 1 )
    {
      //  If interval consist of 2 elements then it's simple
      //
      codes[ptable[li].ch] += '0';
      codes[ptable[ri].ch] += '1';
    }
    else
    {
      //  Calculating sum of probabilities at specified interval
      //
      pfull = 0;
      for( i=li; i<=ri; ++i )
      {
        pfull += ptable[i].p;
      }

      //  Searching center
      //
      p = 0;
      isp = -1; // index of split pos
      phalf = pfull * 0.5f;
      for( i=li; i<=ri; ++i )
      {
        p += ptable[i].p;
        if(p <= phalf) {
          codes[ptable[i].ch] += '0';
        }
        else
        {
          codes[ptable[i].ch] += '1';
          if( isp < 0 ) isp = i;
        }
      }
      
      if( isp < 0 ) isp = li+1;

      //  Next step (recursive)
      //
      EncShannon( li, isp-1 );
      EncShannon( isp, ri );
    }
  }
};

int show_usage() {
  printf("Shannon-Fano coding algorithm"NL);
  printf("by Sergey Tikhonov (st@haqu.net)"NL);
  printf(NL);
  printf("Usage: shannon [OPTIONS] input [output]"NL);
  printf("  The default action is to encode input file."NL);
  printf("  -d\tDecode file."NL);
  printf(NL);
  printf("Examples:"NL);
  printf("  shannon input.txt"NL);
  printf("  shannon input.txt encoded.txt"NL);
  printf("  shannon -d encoded.txt"NL);
  printf(NL);
  exit(0);
}

int main(int argc, char **argv)
{
  int i = 1;
  int dFlag = 0;
  char inputFilename[128];
  char outputFilename[128];

  printf(NL);

  if(i == argc) show_usage();

  if(strcmp(argv[i],"-d") == 0) {
    dFlag = 1;
    i++;
    if(i == argc) show_usage();
  }

  strcpy(inputFilename, argv[i]);
  i++;

  if(i < argc) {
    strcpy(outputFilename, argv[i]);
  } else {
    if(dFlag) {
      strcpy(outputFilename, "decoded.txt");
    } else {
      strcpy(outputFilename, "encoded.txt");
    }
  }

  //  Calling encoding or decoding subroutine
  //
  Coder *coder;
  coder = new Coder;
  assert(coder);
  if(!dFlag) {
    coder->Encode( inputFilename, outputFilename );
  } else {
    coder->Decode( inputFilename, outputFilename );
  }
  delete coder;

  printf(NL);

  return 0;
}
