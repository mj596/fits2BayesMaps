#include <iostream>
#include <CCfits/CCfits>
#include <gsl/gsl_matrix.h>
#include "anaBayes.h"
#include "anaLiMa.h"
#include "anaData.h"
#include "print.h"
#include <boost/progress.hpp>

/* Usage:
   fits2BayesMaps <fits file> configuration number (mode):
   1 - Lima Excess [Excess_LiMa]
   2 - Lima Significance [Sig_LiMa]
   101 - Bayes Mode [Mode_Bayes]
   102 - Bayes Significance [Sig_Bayes]
   104 - Bayes Mode & Significance
   111 - Lima Excess && Bayes Mode
   124 - Lima Significance && Bayes Mode && Bayes Significance
   144 - ALL
   - Lima Excess & Significance
   etc.
*/

class BasicMaps
{
  gsl_matrix *On, *Off, *OnExposure, *OffExposure;
  int size;
  void addMap( std::string type, gsl_matrix* data );
  std::string filename;
  int mode;

public:
  BasicMaps( std::string _filename, int mode=144 );
  ~BasicMaps( ){ }

  void addData( );
  //  void saveData( std::string filename, gsl_matrix* data );
  void saveData( std::string filename, gsl_matrix* data, std::string _extension );
  
  double getOn( int i, int j ) { return gsl_matrix_get( On, i, j ); }
  double getOff( int i, int j ) { return gsl_matrix_get( Off, i, j ); }
  double getOnExposure( int i, int j ) { return gsl_matrix_get( OnExposure, i, j ); }
  double getOffExposure( int i, int j ) { return gsl_matrix_get( OffExposure, i, j ); }

  gsl_matrix* getOn( ) { return On; }
  gsl_matrix* getOff( ) { return Off; }
  gsl_matrix* getOnExposure( ) { return OnExposure; }
  gsl_matrix* getOffExposure( ) { return OffExposure; }

  int getSize( ) { return size; }
};

BasicMaps::BasicMaps( std::string _filename, int mode )
{
  size = 0;
  On = NULL;
  Off = NULL;
  OnExposure = NULL;
  OffExposure = NULL;
  filename = _filename;
  
}

void BasicMaps::addData( )
{
  std::cout << "Opening file: " << filename << std::endl;
  std::auto_ptr<CCfits::FITS> pInfile( new CCfits::FITS( filename,CCfits::Read,true ) );
  std::string type;
  
  if( size == 0 )
    {
      type = "On";
      CCfits::ExtHDU& OnExposure = pInfile->extension( type );
      std::valarray<double> content; 
      OnExposure.read(content);
      size = sqrt( content.size( ) );
      std::cout << "Maps size: " << size << "x" << size << std::endl;
    }

  // On
  type = "On";
  std::cout << "Reading " << type << " map" << std::endl;
  On = gsl_matrix_alloc( size, size );
  CCfits::ExtHDU& OnImage = pInfile->extension( type );
  std::valarray<double> valOn; 
  OnImage.read( valOn );
  for( int i=1;i<=size;i++ )
    for( int j=1;j<=size;j++ )
      gsl_matrix_set( On, i-1, j-1, valOn[(i-1)*size+j] );

  // Off
  type = "Off";
  std::cout << "Reading " << type << " map" << std::endl;
  Off = gsl_matrix_alloc( size, size );
  CCfits::ExtHDU& OffImage = pInfile->extension( type );
  std::valarray<double> valOff; 
  OffImage.read( valOff );
  for( int i=1;i<=size;i++ )
    for( int j=1;j<=size;j++ )
      gsl_matrix_set( Off, i-1, j-1, valOff[(i-1)*size+j] );

  // OnExposure
  type = "OnExposure";
  std::cout << "Reading " << type << " map" << std::endl;
  OnExposure = gsl_matrix_alloc( size, size );
  CCfits::ExtHDU& OnExposureImage = pInfile->extension( type );
  std::valarray<double> valOnExposure; 
  OnExposureImage.read( valOnExposure );
  for( int i=1;i<=size;i++ )
    for( int j=1;j<=size;j++ )
      gsl_matrix_set( OnExposure, i-1, j-1, valOnExposure[(i-1)*size+j] );

  // OffExposure
  type = "OffExposure";
  std::cout << "Reading " << type << " map" << std::endl;
  OffExposure = gsl_matrix_alloc( size, size );
  CCfits::ExtHDU& OffExposureImage = pInfile->extension( type );
  std::valarray<double> valOffExposure; 
  OffExposureImage.read( valOffExposure );
  for( int i=1;i<=size;i++ )
    for( int j=1;j<=size;j++ )
      gsl_matrix_set( OffExposure, i-1, j-1, valOffExposure[(i-1)*size+j] );
}

void BasicMaps::saveData( std::string filename, gsl_matrix* data, std::string _extension )
{
  int lastindex = filename.find_last_of(".");
  string rawname = filename.substr(0, lastindex);
  string filename_save = rawname+"_DerivedMaps.fits";
  std::cout << "Saving " << filename_save << "[" << _extension << "]" << std::endl;
  std::auto_ptr<CCfits::FITS> pFitsWrite(0);
  pFitsWrite.reset( new CCfits::FITS( filename_save, CCfits::Write ) );
  std::vector<long> extAx(2,size);
  CCfits::ExtHDU* imageNew = pFitsWrite->addImage(_extension,DOUBLE_IMG,extAx); 
  
  std::valarray<double> temp(size*size); 

  for( int i=0;i<size;i++ )
    for( int j=0;j<size;j++ )
      temp[i*size+j] = gsl_matrix_get( data, i, j );

  imageNew -> write(1,size*size,temp);
}

int main(int argc, char *argv[])
{
  if( argc == 1 )
    {
      std::cout << "Please specify a filename!" << std::endl;
      exit(0);
    }
  else
    {
      std::string filename = argv[1];
      bool le = 1;
      bool ls = 0;
      bool bm = 0;
      bool bs = 0;
      //      std::cout << "argc " << argc << " argv[2] " << argv[2] << std::endl;
      //      if( std::string(argv[2]) == "lima" )
      //	std::cout << ".... " << std::endl;
      if( argc >= 3 && ( (std::string)argv[2] == "excess" || (std::string)argv[2] == "Excess" || (std::string)argv[2] == "exc" || (std::string)argv[2] == "Exc" ) )
	{ 
	  std::cout << "excess" << std::endl;
	  le = 1;
	  bm = 1;
	}
      else if( argc >= 3 && ( (std::string)argv[2] == "significance" || (std::string)argv[2] == "Significance" || (std::string)argv[2] == "sig" || (std::string)argv[2] == "Sig" ) )
	{
	  std::cout << "sig" << std::endl;
	  ls = 1;
	  bs = 1;
	}
      else if( argc >= 3 && ( (std::string)argv[2] == "lima" || (std::string)argv[2] == "LiMa" || (std::string)argv[2] == "Lima" || (std::string)argv[2] == "liMa" ) )
	{
	  std::cout << "lima" << std::endl;
	  le = 1;
	  ls = 1;
	  std::cout << le << " " << ls << std::endl;
	}
      else if( argc >= 3 && ( (std::string)argv[2] == "Bayes" || (std::string)argv[2] == "bayes" || (std::string)argv[2] == "Bay" || (std::string)argv[2] == "bay" ) )
	{
	  std::cout << "bayes" << std::endl;
	  le = 0;
	  ls = 0;
	  bm = 1;
	  bs = 1;
	}
      else if( argc >= 3 && ( (std::string)argv[2] == "all" || (std::string)argv[2] == "All" || (std::string)argv[2] == "ALL" ) )
	{
	  std::cout << "all" << std::endl;
	  le = 1;
	  bm = 1;
	  ls = 1;
	  bs = 1;
	}
      else if( argc >= 6 )
	{
	  std::cout << "user" << std::endl;
	  le = argv[2];
	  ls = argv[3];
	  bm = argv[4];
	  bs = argv[5];
	}
      
      BasicMaps* data = new BasicMaps( filename );
      data -> addData( );

      if( le || ls ) {
	anaLiMa *alima = new anaLiMa( false );    

	gsl_matrix* lima_excess = gsl_matrix_alloc( data->getSize(), data->getSize() );
	gsl_matrix* lima_significance = gsl_matrix_alloc( data->getSize(), data->getSize() );

	print_info("Running Li&Ma analysis...");
	boost::progress_display show_progress( data->getSize() );
	for( int i=0;i<data->getSize();i++ )
	  {
	    ++show_progress;
	      for( int j=0;j<data->getSize();j++ )
		{
		  alima->set( data->getOn(i,j),
			      data->getOff(i,j),
			      data->getOnExposure(i,j),
			      data->getOffExposure(i,j) );      
		  alima->analyse( );
		  if( le ) { gsl_matrix_set( lima_excess, i, j, alima->getExcess( ) ); }
		  if( ls ) { gsl_matrix_set( lima_significance, i, j, alima->getSignificance( ) ); }
		}
	  }
	
	print_info("Done.");
	
	if( le ) { data -> saveData( filename, lima_excess,"Excess_LiMa" ); }
	if( ls ) { data -> saveData( filename, lima_significance,"Sig_LiMa" ); }
	
	gsl_matrix_free( lima_excess );
	gsl_matrix_free( lima_significance );
	
	delete alima; }
      
      if( bm || bs ) {
	anaBayes *abayes = new anaBayes( false );
	anaLiMa *alima = new anaLiMa( false );    
	
	gsl_matrix* bayes_mode = gsl_matrix_alloc( data->getSize()+100, data->getSize()+100 );
	gsl_matrix* bayes_significance = gsl_matrix_alloc( data->getSize(), data->getSize() );

	double mode = 0.0;
	double bayes_sig = 0.0;
	
	print_info("Running Bayes analysis...");
	boost::progress_display show_progress( data->getSize() );
	for( int i=0;i<data->getSize();i++ )
	  {
	    ++show_progress;
	    for( int j=0;j<data->getSize();j++ )
	      {
		alima->set( data->getOn(i,j),
			    data->getOff(i,j),
			    data->getOnExposure(i,j),
			    data->getOffExposure(i,j) );      
		alima->analyse( );

		abayes->set( data->getOn(i,j),
			     data->getOff(i,j),
			     data->getOnExposure(i,j),
			     data->getOffExposure(i,j) );
		
		if( alima->getExcess() <= 0.0 )
		  {
		    if( bm ) { mode = 0.0; }
		    if( bs ) { bayes_sig = 0.0; }
		  }
		else
		  {
		    if( alima->getSignificance() >= 5.0 )
		      {
			if( bm ) { mode = alima->getExcess(); }
			if( bs ) { bayes_sig = alima->getSignificance(); }
		      }
		    else
		      {
			if( !abayes->analyse( false ) )
			  {
			    if( bm ) { mode = abayes->getMode(); }
			    if( bs ) { bayes_sig = abayes->getSignificance(); }
			  }
			else
			  {
			    if( bm ) { mode = 0.0; }
			    if( bs ) { bayes_sig = 0.0; }
			  }
		      }
		  }
		if( bm ) { gsl_matrix_set( bayes_mode, i, j, mode ); }
		if( bs ) { gsl_matrix_set( bayes_significance, i, j, bayes_sig ); }
	      }
	  }
	
	print_info("Done.");
	
	if( bm ) { data -> saveData( filename, bayes_mode, "Mode_Bayes" ); }
	if( bs ) { data -> saveData( filename, bayes_significance, "Sig_Bayes" ); }
	
	gsl_matrix_free( bayes_mode );
	gsl_matrix_free( bayes_significance );
      
	delete abayes;
	delete alima;
      }
    }
  return 0;
}

