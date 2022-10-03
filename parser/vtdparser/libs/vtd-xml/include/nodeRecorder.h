/* 
 * Copyright (C) 2002-2015 XimpleWare, info@ximpleware.com
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
/*VTD-XML is protected by US patent 7133857, 7260652, an 7761459*/
#pragma once
#ifndef NR_H
#define NR_H
#include "vtdNav.h"
namespace vtdxml{
	/**
	 * @brief Node record allows one to record the node position of VTDNav. 
	 * You can store/load multiple node position using NodeRecorder, 
	 * which is also more spatially efficient than VTDNave. However, the internal 
	 * representation of node is variable in length and recording a large 
	 * number of nodes could consume much memory. It is a bad idea to record
	 *  every node of an document. So be careful when using it.
	 */
	class NodeRecorder{
	public:
		NodeRecorder();
		NodeRecorder(VTDNav *vn);
		virtual ~NodeRecorder();

		/**
		 * @brief Record the position of the VTDNav into the internal buffer
		 * 
		 */
		void record();

		/**
		 * @brief Set the pointer to the first node in NodeRecorder.
		 * 
		 */
		void resetPointer();

		/**
		 * @brief Erase all the recorded nodes in NodeRecorder and
		 * recycle the internal buffer.
		 * 
		 */
		void clear();

		/**
		 * @brief This method set the cursor in VTDNav to the nodes as recorded
		 * in NodeRecorder, and return the output of "getCurrentIndex()".
		 * It is important to notice that you can only go forward, not backward.
		 * 
		 * @return int 
		 */
		int iterate();

		/**
		 * @brief Bind the NodeRecorder to a VTDNav object.
		 * 
		 * @param vn 
		 */
		void bind(VTDNav *vn);
		static const int BUF_SZ_EXPO  = 7;
	
	private:
		VTDNav *vn;
		// internal buffer for storing the pointers
		FastIntBuffer *fib;
		int position;
		int size;
		int count;
	};
}
#endif