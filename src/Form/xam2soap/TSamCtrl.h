/*
 ******************************************************************************
 *Copyright 2010
 *	BGI-SHENZHEN
 *All Rights Reserved
 *ATHUOR : Bill Tang
 *CREATE DATE : 2010-7-29
 *CLASS NAME: TSamCtrl
 *FUNCTION : a class that can control sam/bam format file. It is convenient to control sam/bam 
 *			format file like open, read, write, close and so on.
 *FILE NAME : TSamCtrl.h
 *UPDATE DATE : 2010-7-30
 *UPDATE BY : Bill Tang
 *******************************************************************************
 */

#ifndef TTSamCtrl_H_
#define TTSamCtrl_H_ 


#pragma once

#include <cmath>
#include <string>
#include "Ttools.h"
#include  <zlib.h>
#include <htslib/sam.h>
//#include <sam_view.h>

//#include "../../include/sam/bam.h"

class TSamCtrl
{
	private:
		std::string m_in_path; // the sam/bam in file's path.
		std::string m_out_path; // the sam/bam out file's path.
		samfile_t *m_in; // sam/bam file handler, use as reading.
		samfile_t *m_out; // sam/bam file handler, use as writing.
		std::string m_in_mode; // reading mode
		std::string m_out_mode; // wrinting mode
		char *m_fn_list; // sam/bam file list.
		char *m_s; // an tmp varible
		bam1_t *m_b; // a bam constructor.
	public:
		TSamCtrl();
		~TSamCtrl(void);

		bool open(const char *path, const char *mode);
		void close();
		int readline(std::string &line);
		bool isOpened();
};


/*
 *****************************************************************************
 *Copyright 2010
 *	BGI-SHENZHEN
 *All Rights Reserved
 *ATHUOR : Bill Tang
 *CREATE DATE : 2010-7-29
 *FUNCTION : definition of class TSamCtrl
 *FILE NAME : TSamCtrl.cpp
 *UPDATE DATE : 2010-8-1
 *UPDATE BY : Bill Tang
 *******************************************************************************
 */


//#include "TSamCtrl.h"


TSamCtrl::TSamCtrl()
{
	m_out_path = "-";
	m_out_mode = "w";
	m_fn_list = 0;
	m_in = 0;
	m_out = 0;
	m_b = bam_init1(); // initialize m_b

}


TSamCtrl::~TSamCtrl(void)
{
	bam_destroy1(m_b);
	close();
}

/**
 * DATE: 2010-7-29
 * FUNCTION: open sam/bam file in pattern mode.
 * PARAMETER: path: the file's path. mode: the pattern.
 * RETURN: successful: true. failed: flase.
 */
bool TSamCtrl::open(const char *path, const char *mode) 
{
	if (mode[0] == 'r') {
		// judge if opened a file before
		if (m_in != 0) {
			return false;
		}

		m_in_mode = mode;
		m_in_path = path;

		// open sam/bam file
		if ((m_in = samopen(m_in_path.c_str(), m_in_mode.c_str(), m_fn_list)) == 0) 
		{
			return false;		
		}

		if (m_in->header == 0) 
		{
			return false;
		}

		if ((m_out = samopen(m_out_path.c_str(), m_out_mode.c_str(), m_in->header)) == 0) 
		{
			return false;
		}		
		return true;
	}
	else
	{
		return false;
	}
}

/**
 * DATE: 2010-7-29
 * FUNCTION: read a line from sam/bam file.
 * PARAMETER: line: the read data will stored in this varible.
 * RETURN: the line size when read successful. -1 when read failed
 */
int TSamCtrl::readline(std::string &line) 
{
	if (m_in == 0) 
	{
		if (m_fn_list != 0)
			delete m_fn_list;
		if (m_in != 0)
			samclose(m_in);
		if (m_out != 0)
			samclose(m_out);	    

		m_fn_list = 0;
		m_out = 0;
		return -1;
	}

	int ret = 0;
	// begin to read
	while ((ret = samread(m_in, m_b)) >= 0) 
	{
		// when read failed continue
		if (__g_skip_aln(m_in->header, m_b)) 
		{
			continue;
		}

		m_s = bam_format1_core(m_out->header, m_b, m_out->type>>2&3); // read the buffer
		line = m_s; // store into the line
		free(m_s);
		return line.size();
	}

	//if (m_fn_list != 0){  delete m_fn_list;}
	//if (m_in != 0){ samclose(m_in);}
	//if (m_out != 0) { samclose(m_out); }
	m_in = 0;
	m_out = 0;
	m_fn_list = 0;

	return -1; 
}

/**
 * DATE: 2010-7-29
 * FUNCTION: close the sam/bam file
 * PARAMETER: void
 * RETURN: void
 */
void TSamCtrl::close() 
{
	if (m_fn_list != 0)
		delete m_fn_list;
	if (m_in != 0)
		samclose(m_in);
	if (m_out != 0)
		samclose(m_out);
	m_in = 0;
	m_out = 0;
	m_fn_list = 0;
}

/**
 * DATE: 2010-7-29
 * FUNCTION: judge if the file is opened successful
 * PARAMETER: void
 * RETURN: true if successful else false if failed
 */
bool TSamCtrl::isOpened() 
{
	return !(m_in == 0);
}

#endif   //


