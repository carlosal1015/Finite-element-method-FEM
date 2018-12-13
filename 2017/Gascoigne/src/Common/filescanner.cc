/**
*
* Copyright (C) 2004, 2005, 2007, 2011 by the Gascoigne 3D authors
*
* This file is part of Gascoigne 3D
*
* Gascoigne 3D is free software: you can redistribute it
* and/or modify it under the terms of the GNU General Public
* License as published by the Free Software Foundation, either
* version 3 of the License, or (at your option) any later
* version.
*
* Gascoigne 3D is distributed in the hope that it will be
* useful, but WITHOUT ANY WARRANTY; without even the implied
* warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
* PURPOSE.  See the GNU General Public License for more
* details.
*
* Please refer to the file LICENSE.TXT for further information
* on this license.
*
**/


#include  "filescanner.h"
#include  "fixarray.h"
#include  "linescanner.h"
#include  "stlio.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>

using namespace std;

namespace Gascoigne
{
  
/***************************************************/

FileScanner::FileScanner(DataFormatHandler& D) : DH(D)
{
  complain = 1;
  blocksymbol                        = "//Block";
  _i_defaultvalues_level             = 0;
  _i_defaultvalues_save_all_to_file  = 0;
  _s_defaultvalues_save_filename     = "allparamfilevalues.pl";
}

/***************************************************/

FileScanner::FileScanner(DataFormatHandler& D, const ParamFile* pf, const string& blockname) : DH(D)
{
  complain    = 1;
  blocksymbol = "//Block";
  _i_defaultvalues_level             = 0;
  _i_defaultvalues_save_all_to_file  = 0;
  _s_defaultvalues_save_filename     = "allparamfilevalues.pl";
  readfile(pf,blockname);
}

/***************************************************/

void FileScanner::_assert(bool b, const vector<string>& words) const
{
  if(!b) {
    cerr << "*** FileScanner:\tWrong number of arguments in row\t" << words << endl;
    abort();
  }
}

/***************************************************/

void FileScanner::readfile(const ParamFile* pf, const string& blockname)
{
  if(pf==NULL)
    {
      return;
    }


  string inputname = pf->GetName();
  LineScanner LS(inputname);

  if(blockname=="Multilevelsolver"){
    std::cerr << "Error: Trying to search file '<< inputname <<' for Blockname 'Multilevelsolver'."<<std::endl;
    std::cerr << "The blockname 'Multilevelsolver' is not used any longer! It has changed to MultiLevelSolver."<<std::endl;
    std::cerr << "Update code and file."<<std::endl;
    abort();
  }

  vector<string> words;
  int nwords = 0;
  bool searchblock = (blockname!="");
  bool blockfound  = !searchblock;
  bool helpfound   = 0;
  string helpname = "HELP_ME";

  // der folgende block macht es moeglich beliebig viele defaultwert-dateien anzugeben
  // als massnahme gegen schleifen der inklusion gibt es den defaultvalues-zaehler
  // bei jedem weiteren auslesen der defaultwerte der 'naechst-tiefen' datei wird der zaehler
  // um eins hoch gezaehlt
  //syntax: 
  //   //Block DefaultValues
  //   files   3   file_a  file_b  file_c
  //   file_a  settingsa.param
  //   file_b  settingsb.param
  //   file_c  settingsc.param
  //
  if(blockname!="DefaultValues"){
     _i_defaultvalues_level++;
     if( 10 < _i_defaultvalues_level ){
       cerr<< __FILE__ << ":" << __LINE__ << ": there seems to be a 'Block DefaultValues' loop"<<endl;
       cerr<< __FILE__ << ":" << __LINE__ << ": last used file: "<< pf->GetName() <<endl;
       abort();
     }
     DataFormatHandler DFH;    
     vector<string>    vs_files;
     string            s_paramfile = pf->GetName();

     DFH.insert("files"            , &vs_files                                                              );
     DFH.insert("save_all_to_file" , &_i_defaultvalues_save_all_to_file , _i_defaultvalues_save_all_to_file );
     DFH.insert("save_filename"    , &_s_defaultvalues_save_filename    , _s_defaultvalues_save_filename    );
     FileScanner FS(DFH);
     FS._i_defaultvalues_level = _i_defaultvalues_level;
     FS.NoComplain();
     FS.readfile(pf,"DefaultValues");

     for(int i=0;i<vs_files.size();i++){
       DataFormatHandler DFH2;
       string            s_filename;
       string            s_keyname = vs_files[i];

       DFH2.insert(s_keyname,&s_filename,"none");  
       FileScanner FS(DFH2);
       FS._i_defaultvalues_level = _i_defaultvalues_level;
       FS.NoComplain();
       FS.readfile(pf,"DefaultValues");

       if(s_filename!="none" && s_filename!= s_paramfile){
         ParamFile paramfile(s_filename);
         FileScanner FS2(DH);
         FS2._i_defaultvalues_level = _i_defaultvalues_level;
         FS2.NoComplain();
         FS2.readfile(&paramfile,blockname);
       }
     }
  }

  while (nwords>=0)
    {
      nwords = LS.NextLine(words);

      if (nwords==0) continue;

      if (words[0]!=blocksymbol) continue;

      if (nwords==1)
        {
          cout << "FileScanner::Block without name" << endl;
        }
      else if (words[1]==helpname)
        {
          helpfound = 1;
        }
      else if (words[1]==blockname)
        {
          blockfound = 1;
          break;
        }
    }
  if (helpfound) print(blockname);

  // this has to go if the save functionality at the bottom is to work
  // if (!blockfound)
  //   {
  //     //cout << "FileScanner::missing Block " << blockname << endl;
  //     return;
  //   }
  //
  // scanning parameters in block
  //
  if (blockfound) {
    while (nwords>=0) {
      nwords = LS.NextLine(words);

      if (nwords==0) continue;

      // testing if next block begins
      //
      if (words[0]==blocksymbol) break;
      //
      // testing commentaries
      //
      if (words[0]=="/*") continue;
      if (words[0]=="//") continue;

      if (nwords==1) {
        cout << "where is the parameter \"" << words[0] << "\" ?" << endl;
        continue;
      }
      FormatToValue(words);
    }
  }

  if( _i_defaultvalues_level==1  && _i_defaultvalues_save_all_to_file){
    // By setting the the paramfile values
    //
    //   //Block DefaultValues
    //   save_all_to_file 1
    //   save_filename    file.pl  // (default is allparamfilevalues.pl)
    //
    // one can save *all* paramfile values that are used by Gascoigne. This includes
    // also the default values of options that the user did not explicitly set in the paramfile.
    // (It also includes values that may have been loaded in other include files. The
    // DefaultValues Block is not saved, since the include commands don't belong in a 'flattened file')
    // This is a good and easy way of determining and listing all availible
    // Gascoigne options.
    //
    // The generated file is a perl script that has to be executed to
    // create the paramfile.

    fstream fcheckforfile(_s_defaultvalues_save_filename.c_str(),ios::in);
    if( !fcheckforfile.is_open() ) {
      // if the file does not yet exist -> create perl header of the file
      // the perl script reads "itself" and generates the complete paramfile automatically (it writes to stdout);
      // the values and blocks may occur multiply in the file, only the last entries are used
      ofstream createheader(_s_defaultvalues_save_filename.c_str() , ios_base::out);
      createheader << "#!/usr/bin/env perl"                                           <<endl;
      createheader << "# to generate param-file, simply execute this perl-script"     <<endl;
      createheader << "$colwidth=35;"                                                 <<endl;
      createheader << ""                                                              <<endl;
      createheader << "open(FH, $0);"                                                 <<endl;
      createheader << "$blockname='nix';"                                             <<endl;
      createheader << "while(<FH>){"                                                  <<endl;
      createheader << "  if(m|^<Block\\s+name='([^']+)'\\s+|){  # '"                  <<endl;
      createheader << "    $blockname=$1; last;"                                      <<endl;
      createheader << "  }"                                                           <<endl;
      createheader << "}"                                                             <<endl;
      createheader << "%data=(); "                                                    <<endl;
      createheader << "while(<FH>){"                                                  <<endl;
      createheader << "  if(m|^<Block\\s+name='([^']+)'\\s+|){ # '"                   <<endl;
      createheader << "    $blockname=$1;"                                            <<endl;
      createheader << "  }elsif(m|^</Block>|){"                                       <<endl;
      createheader << "    $blockname='nix';"                                         <<endl;
      createheader << "  }elsif(m|^\\s*([^\\s]+)\\s+([^\\s].*)\\s*$|){"               <<endl;
      createheader << "    $data{\"$blockname:PERLSPLIT:$1\"}=$2;"                    <<endl;
      createheader << "  }"                                                           <<endl;
      createheader << "}"                                                             <<endl;
      createheader << "$blockname='nix';"                                             <<endl;
      createheader << "foreach $key (sort(keys %data )){"                             <<endl;
      createheader << "  (@names) = split(':PERLSPLIT:',$key);"                       <<endl;
      createheader << "  if($names[0] ne $blockname){"                                <<endl;
      createheader << "    $blockname= $names[0];"                                    <<endl;
      createheader << "    print \"\\n//Block $blockname\\n\";"                       <<endl;
      createheader << "  }"                                                           <<endl;
      createheader << "  if($names[1]=~m|^(.+)\\[(.+)\\]$|){"                         <<endl;
      createheader << "    $names[1]=$1;"                                             <<endl;
      createheader << "    $data{$key}=\"$2     $data{$key}\";"                       <<endl;
      createheader << "  }"                                                           <<endl;
      createheader << "  printf(\"\\%-${colwidth}s \\%s\\n\",$names[1],$data{$key});" <<endl;
      createheader << "}"                                                             <<endl;
      createheader << "print \"\\n//Block nix\\n\";"                                  <<endl;
      createheader << "close(F);"                                                     <<endl;
      createheader << "__END__"                                                       <<endl;
      createheader.close();

      chmod(_s_defaultvalues_save_filename.c_str(), 00755); 
    }else{
      fcheckforfile.close();
    }

    // with this the date is added as a tag for each saved block;
    // this way the file can be used to achive the history of changes in a paramfile;
    // the perl script only uses the last changes in the file; the date tag is
    // ignored by the script
    char ca_date[100];
    {
      time_t  t_date;
      struct tm *tmzgr;
      t_date = time(NULL);
      tmzgr = localtime(&t_date);
      strftime(ca_date,100,"%Y.%m.%d-%H:%M",tmzgr); 
    }

    ofstream savefile(_s_defaultvalues_save_filename.c_str() , ios_base::app);
    savefile << "<Block name='"<< blockname << "' file='"<< pf->GetName() <<"' date='"<< ca_date <<"'>\n";
    DH.print(savefile);
    savefile << "</Block>\n\n";
    savefile.close();
  }
}

/***************************************************/

void FileScanner::FormatToValue(const vector<string>& words)
{
  string keyword_type;
  string keyword = words[0];
  DH.get(keyword_type,keyword);

  /*----------------------------------------------*/

  if (keyword_type=="string")
    {
      DH.setvalue(keyword,words[1]);
    }
  else if (keyword_type=="integer")
    {
      int value = atoi(words[1].c_str());
      DH.setvalue(keyword,value);
    }
  else if (keyword_type=="bool")
    {
      if(words[1]=="true" || words[1]=="1")
        {
          DH.setvalue(keyword,true);
        }
      else
        {
          DH.setvalue(keyword,false);
        }
    }
  else if (keyword_type=="float")
    {
      float value = atof(words[1].c_str());
      DH.setvalue(keyword,value);
    }
  else if (keyword_type=="double")
    {
      double value = atof(words[1].c_str());
      DH.setvalue(keyword,value);
    }

  /*----------------------------------------------*/

  else if (keyword_type=="fixarray<2,double>")
    {
      fixarray<2,double> value;
      value[0] = atof(words[1].c_str());
      value[1] = atof(words[2].c_str());
      DH.setvalue(keyword,value);
    }

  /*----------------------------------------------*/

  else if (keyword_type=="vector<double>")
    {
      int n = atoi(words[1].c_str());
      _assert(words.size()>n+1,words);
      // wenn n==0, dann ist auch n=0 gemeint: rufe setvalue auch dafuer auf
      if (n>=0)
        {
          vector<double> value(n);
          for (int i=0; i<n; i++)
            {
              value[i] = atof(words[i+2].c_str());
            }
          DH.setvalue(keyword,value);
        }
    }
  else if (keyword_type=="IntVector")
    {
      int n = atoi(words[1].c_str());
      _assert(words.size()>n+1,words);
      // cerr << "-----" << n << " " << words.size() << endl;
      // wenn n==0, dann ist auch n=0 gemeint: rufe setvalue auch dafuer auf
      if (n>=0)
        {
          IntVector value(n);
          for (int i=0; i<n; i++)
            {
              value[i] = atoi(words[i+2].c_str());
            }
          DH.setvalue(keyword,value);
        }
    }
  else if (keyword_type=="vector<string>")
    {
      int n = atoi(words[1].c_str());
      _assert(words.size()>n+1,words);
      // wenn n==0, dann ist auch n=0 gemeint: rufe setvalue auch dafuer auf
      if (n>=0)
        {
          vector<string> value(n);
          for (int i=0; i<n; i++)
            {
              value[i] = words[i+2];
            }
          DH.setvalue(keyword,value);
        }
    }

  /*----------------------------------------------*/

  else if (keyword_type=="map<int,IntVector >")
    {
      int col = atoi(words[1].c_str());
      int n   = atoi(words[2].c_str());
      IntVector value(n);
      _assert(words.size()>n+2,words);
      for(int i=0; i<n; i++) value[i] = atoi(words[i+3].c_str());
      pair<int,IntVector > p = make_pair(col,value);
      DH.setvalue(keyword,p);
    }

  /*----------------------------------------------*/

  else if (keyword_type=="set<vector<string> >")
    {
      int n = atoi(words[1].c_str());
      vector<string> value(n);
      _assert(words.size()>n+1,words);
      for(int i=0; i<n; i++) value[i] = words[i+2];

      DH.insertvalue(keyword,value);
    }

  else if (keyword_type=="set<int>")
    {
      int n = atoi(words[1].c_str());
      set<int> value;
      _assert(words.size()>n+1,words);
      for(int i=0; i<n; i++) 
        {
          int v = atoi(words[i+2].c_str());
          value.insert(v);
        }

      DH.setvalue(keyword,value);
    }

  /*----------------------------------------------*/

  else if (keyword_type=="StringDouble")
    {
      pair<string,double> p = make_pair(words[1], atof(words[2].c_str()));
      DH.setvalue(keyword,p);
    }
  else if (keyword_type=="" && complain)
    {
      cerr << " ********  FileScanner:: keyword '"<< keyword <<"' found in param-file but not registered in the DataFormatHandler-code."<< endl;
    }
}

/***************************************************/

void FileScanner::print(const string& blockname) const
{
  cout << "=====================" << endl;
  cout << blocksymbol << " " << blockname << endl << endl;
  DH.print(cout);
}
}
