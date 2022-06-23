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
#ifndef ELEMENT_FRAGMENT_NS_H
#define ELEMENT_FRAGMENT_NS_H
//#include "customTypes.h"
#include "vtdNav.h"
//#include "fastIntBuffer.h"

namespace com_ximpleware {
	/**
	 * @brief A namespace compensated element fragment.
	 */
	class ElementFragmentNs{
	private:
		Long l;
		FastIntBuffer *fib;
		int stLen;
		VTDNav* vn;
		//bool ws_init = false;
        
	public:
		ElementFragmentNs(VTDNav *vn, Long l1, FastIntBuffer *fib, int len);
		virtual ~ElementFragmentNs();

		/**
		 * @brief Get the byte length of the namespace compensated element fragment
		 * in its source encoding format.
		 * 
		 * @return int 
		 */
		int getFragmentSize();

		/**
		 * @brief Get the byte length of the namespace compensated element fragment
		 * in the specified encoding format.
		 * 
		 * @param t 
		 * @return int the target encoding
		 */
		int getFragmentSize(encoding_t t);

		/**
		 * @brief Get the long encoding of the len and offset of the element fragment.
		 * 
		 * @return Long the long encoding
		 */
		Long getOffsetLeng();
		
		/**
		 * @brief Convert to a byte array of an element with namespace compensation in
		 * the source encoding format.
		 * 
		 * @return UByte* the byte array
		 */
		UByte* toFragmentBytes();

		/**
		 * @brief Convert to a byte array of an element with namespace compensation in
		 * the specified encoding format.
		 * 
		 * @param dest_encoding the target encoding
		 * @return UByte* the byte array
		 */
		UByte* toFragmentBytes(encoding_t dest_encoding);

		/**
		 * @brief Write the namespace compensated element fragment to the FILE pointer.
		 * 
		 * @param f the FILE pointer to write to
		 */
		void writeFragmentToFile(FILE *f);

		/**
		 * @brief Write the namespace compensated element fragment to the FILE pointer
		 * in the specified encoding format.
		 * 
		 * @param f the FILE pointer to write to
		 * @param t the target encoding
		 */
		void writeFragmentToFile(FILE *f, encoding_t dest_encoding);

		static UByte ws[5];
	};
};
#endif