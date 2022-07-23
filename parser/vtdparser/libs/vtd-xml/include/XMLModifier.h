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
#ifndef XMLMODIFIER_H
#define XMLMODIFIER_H
#include "vtdGen.h"
#include "elementFragmentNs.h"
#include "intHash.h"
#include "transcoder.h"

namespace vtdxml {
	/**
	 * @brief XMLModifier offers an easy-to-use interface for users to 
	 * take advantage of the incremental update of VTD-XML.
	 * The XML modifier assumes there is a master document on which the 
	 * modification is applied: users can remove an element, update a token, 
	 * replace an element name, or insert new content anywhere in the document
	 * transcoding methods are built-in.
	 * 
	 * Note: The modification operations are recorded first then, output() is called 
	 * to generate output document
	 */
	class XMLModifier {
		
	public:
		const static Long MASK_DELETE = 0x00000000000000000LL;
		const static Long MASK_INSERT_SEGMENT_BYTE  =0x2000000000000000LL;
		const static Long MASK_INSERT_BYTE =0x4000000000000000LL;
		const static Long MASK_INSERT_SEGMENT_STRING =0x6000000000000000LL;
		const static Long MASK_INSERT_STRING =0x8000000000000000LL;
		const static Long MASK_INSERT_FRAGMENT_NS  =0xa000000000000000LL;
		const static Long MASK_INSERT_SEGMENT_BYTE_ENCLOSED = 0x6000000000000000LL;
		const static Long MASK_INSERT_BYTE_ENCLOSED = 0x8000000000000000LL;
		const static Long MASK_INSERT_FRAGMENT_NS_ENCLOSED = 0xe000000000000000LL;
		const static Long MASK_NULL =  0xc000000000000000LL; //1100

		XMLModifier();
		XMLModifier(VTDNav *vn);
		~XMLModifier();

		/**
		 * @brief Attach master document to this XMLModifier so that all further operations
		 * are performed on the master document.
		 * 
		 * @param md master document
		 */
		void bind(VTDNav *md);

		/**
		 * @brief Remove content from the master document. It first calls getCurrentIndex();
		 * if the result is a starting tag, then the entire element referred to by the starting tag
		 * is removed; if the result is an attribute name or ns node, then the corresponding attribute
		 * or value pair is removed. If the token type is text, CDATA, or comment, then the entire node
		 * including the starting and ending delimiting text is removed.
		 * 
		 */
		void remove();

		/**
		 * @brief Remove the token content at index i
		 * 
		 * @param i 
		 */
		void removeToken(int i);

		/**
		 * @brief Remove the attribute value pair at index i
		 * 
		 * @param attrNameIndex 
		 */
		void removeAttribute(int attrNameIndex);

		/**
		 * @brief Remove a segment of byte content of size len
		 * 
		 * @param offset 
		 * @param len 
		 */
		void removeContent(int offset, int len);

		/**
		 * @brief Update the token with the given Unicode string.
		 * 
		 * @param index 
		 * @param newContent 
		 */
		void updateToken( int index, UCSChar *newContent);

		/**
		 * @brief Update the token with the given byte array at the given position with
		 * offset.
		 * 
		 * @param index 
		 * @param newContent 
		 */
		void updateToken(int index, UByte *byteContent, int contentOffset, int contentLen);

		/**
		 * @brief Update the token with the transcoded byte array at the given position with
		 * offset.
		 * 
		 * @param index 
		 * @param byteContent 
		 * @param contentOffset 
		 * @param contentLen 
		 * @param src_encoding 
		 */
		void updateToken(int index, UByte *byteContent, int contentOffset, int contentLen, encoding_t src_encoding);

		/**
		 * @brief Update the token with the transcode representation of
		 * a segment of a byte array contained in the VTDNav instanct vn
		 * at the given position with offset.
		 * 
		 * @param index 
		 * @param vn
		 * @param contentOffset
		 * @param contentLen
		 */
		void updateToken(int index, VTDNav *vn, int contentOffset, int contentLen);

		/**
		 * @brief Insert the Unicode string after the element.
		 * 
		 * @param s a UCSChar string
		 */
		void insertAfterElement(UCSChar *s);

		/**
		 * @brief Insert the Unicode string before the element.
		 * 
		 * @param s a UCSChar string
		 */
		void insertBeforeElement(UCSChar *s);

		/**
		 * @brief Insert the Unicode string as an attribute after the starting tag.
		 * 
		 * @param s a UCSChar string
		 */
		void insertAttribute(UCSChar *attr);

		/**
		 * @brief Insert the Unicode string after the head of the cursor.
		 * 
		 * @param attr 
		 */
		void insertAfterHead(UCSChar *s);

		/**
		 * @brief This method will first call getCurrentIndex() to get the cursor index value
		 * then insert the byte array ba after the element
		 * 
		 * @param ba the byte array to be inserted
		 * @param arrayLen length of the byte array
		 */
		void insertAfterElement(UByte* ba, int arrayLen);

		/**
		 * @brief This method will first call getCurrentIndex() to get the cursor index value
		 * then insert the character array ba after the element
		 * 
		 * @param ba the char array to be inserted
		 * @param arrayLen length of the char array
		 */
		void insertAfterElement(char* ba,int arrayLen){
			insertAfterElement((UByte *)ba,arrayLen);
		}

		void insertBeforeElement(UByte* ba, int arrayLen);
		void insertBeforeElement(char* ba,int arrayLen){
			insertBeforeElement((UByte *)ba,arrayLen);
		}
		void insertAfterHead(UByte* ba, int arrayLen);
		void insertAfterHead(char* ba,int arrayLen){
			insertAfterHead((UByte *)ba,arrayLen);
		}

		void insertAfterElement(UByte* ba, int contentOffset, int contentLen);
		void insertAfterElement(char* ba, int contentOffset, int contentLen){
			insertAfterElement((UByte *)ba, contentOffset,contentLen);
		}
		void insertBeforeElement(UByte* ba, int contentOffset, int contentLen);
		void insertBeforeElement(char* ba, int contentOffset, int contentLen){
			insertBeforeElement((UByte *)ba, contentOffset,contentLen);
		}
		void insertAfterHead(UByte* ba, int contentOffset, int contentLen);
		void insertAfterHead(char* ba, int contentOffset, int contentLen){
			insertAfterHead((UByte *)ba, contentOffset,contentLen);
		}
		void insertBeforeElement(ElementFragmentNs *ef);
		void insertAfterElement(ElementFragmentNs *ef);
		void insertAfterHead(ElementFragmentNs *ef);


		void insertAfterElement(encoding_t src_encoding, UByte* ba, int arrayLen);
		void insertAfterElement(encoding_t src_encoding, char* ba, int arrayLen){
			insertAfterElement(src_encoding, (UByte *) ba, arrayLen);
		}
		void insertBeforeElement(encoding_t src_encoding, UByte* ba, int arrayLen);
		void insertBeforeElement(encoding_t src_encoding, char* ba, int arrayLen){
			insertBeforeElement(src_encoding, (UByte *) ba, arrayLen);
		}
		void insertAfterHead(encoding_t src_encoding, UByte* ba, int arrayLen);
		void insertAfterHead(encoding_t src_encoding, char* ba, int arrayLen){
			insertAfterHead(src_encoding, (UByte *) ba, arrayLen);
		}

		void insertAfterElement(encoding_t src_encoding, UByte* ba, int contentOffset, int contentLen);
		void insertBeforeElement(encoding_t src_encoding, UByte* ba, int contentOffset, int contentLen);
		void insertAfterHead(encoding_t src_encoding, UByte* ba, int contentOffset, int contentLen);

		
		void insertAfterElement(encoding_t src_encoding, char* ba, int contentOffset, int contentLen){
			insertAfterElement(  src_encoding, (UByte*) ba, contentOffset, contentLen);
		}
		void insertBeforeElement(encoding_t src_encoding, char* ba, int contentOffset, int contentLen){
			insertBeforeElement( src_encoding, (UByte*) ba, contentOffset, contentLen);
		}
		void insertAfterHead(encoding_t src_encoding, char* ba, int contentOffset, int contentLen){
			insertAfterHead(  src_encoding, (UByte*) ba, contentOffset, contentLen);
		}


		void insertAfterElement(VTDNav *vn1, int contentOffset, int contentLen);
		void insertBeforeElement(VTDNav *vn1, int contentOffset, int contentLen);
		void insertAfterHead(VTDNav *vn1, int contentOffset, int contentLen);

		///2.11----------------------
		void insertBeforeTail(VTDNav *vn1, int contentOffset, int contentLen){
			insertBeforeTail(vn1->encoding,vn1->XMLDoc,contentOffset, contentLen); 
		}
		/*void insertBeforeTail(VTDNav *vn1, Long l){
			insertBeforeTail(vn1->encoding, vn1->XMLDoc, l);
		}*/
		
		/*void insertBeforeTail(UByte *ba, Long l);
		void insertBeforeTail(char *ba, Long l){
			insertBeforeTail((UByte *)ba, l);
		}*/
		
		void insertBeforeTail(UByte *ba, int contentOffset, int contentLen);
		void insertBeforeTail(char *ba, int contentOffset, int contentLen){
			insertBeforeTail((UByte *)ba, contentOffset, contentLen);
		}

		void insertBeforeTail(UByte *ba, int arrayLen);
		void insertBeforeTail(char *ba, int arrayLen){
			insertBeforeTail((UByte *)ba, arrayLen);
		}

		void insertBeforeTail(ElementFragmentNs *ef);

		/*void insertBeforeTail(encoding_t src_encoding,UByte *ba, Long l);
		void insertBeforeTail(encoding_t src_encoding,char *ba, Long l){
			insertBeforeTail(src_encoding, (UByte *)ba, l);
		}*/

		void insertBeforeTail(encoding_t src_encoding, UByte *ba, int contentOffset, int contentLen);
		void insertBeforeTail(encoding_t src_encoding,char *ba,int contentOffset, int contentLen){
			insertBeforeTail(src_encoding, (UByte *)ba, contentOffset, contentLen);
		}

		void insertBeforeTail(encoding_t src_encoding, UByte *ba, int arrayLen);
		void insertBeforeTail(encoding_t src_encoding, char *ba, int arrayLen){
			insertBeforeTail(src_encoding, (UByte *)ba, arrayLen);
		}

		void insertBeforeTail(UCSChar *s);

		/*
		void updateToken2(XMLModifier *xm, int index, UByte *newContentBytes, int len);
		void insertAfterElement2(XMLModifier *xm, UByte *b, int len);
		void insertBeforeElement2(XMLModifier *xm, UByte *b, int len);
		void insertAttribute2(XMLModifier *xm, UByte *attr, int len);
		*/

		void output(FILE *f);
		void output(char *fileName);

		void updateElementName(UCSChar* elementName);

		void reset();

	private:
		
		typedef Long (XMLModifier::*getBytes)(UCSChar *s);
		encoding_t encoding;
		IntHash *deleteHash;
		IntHash *insertHash;
		FastLongBuffer *flb; /* lower 32 bit offset, upper 29 bits
						 length, upper 3 bits */
		FastLongBuffer *fob; /*lower 32 bits the object pointer,
						 upper 32 bits the length of the byte array*/
		VTDNav *md; /*master document*/
		getBytes gbytes;

		Long getBytes_UTF8(UCSChar *s);
		Long getBytes_UTF16LE(UCSChar *s);
		Long getBytes_UTF16BE(UCSChar *s);
		Long getBytes_ISO_8859_1(UCSChar *s);
		Long getBytes_ASCII(UCSChar *s);

		void check( );
		void sort( );
		void quickSort( int lo, int hi);
		void insertBytesAt( int offset, Long l);
		void insertBytesAt2( int offset, Long lenPlusPointer);
		void insertBytesAt( int offset, ElementFragmentNs* ef);
		UByte *doubleCapacity(UByte *b, size_t cap);

		void insertBytesEnclosedAt(int offset, Long l);
		void insertBytesEnclosedAt2(int offset, Long lenPlusPointer);
		void insertBytesEnclosedAt( int offset, ElementFragmentNs* ef);

		void insertEndingTag(Long l);
		void check2();
		
		
	};
}

#endif