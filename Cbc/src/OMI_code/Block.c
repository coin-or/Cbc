/**
  * @(#)Block.c
  * @author Robert Entriken<br>
  * Copyright (C) 1999-2007 EPRI<br>
  * All rights reserved.
  * @version 07-03-07
  * @since OMI_1.0
  *
  * Revisions:
  *   02-06-16 bobe     Switched to MW generated headers
  *   06-02-20 bobe     revised copyright
  *   06-01-28 jmm      new package layout + casts/corrections
  *   06-04-10 jmm      cast NULL to jint and some jint to integer to avoid warnings
  *   07-03-07 bobe     updated Copyright
  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "omi_Block.h"
#include "com_epri_omi_Base.h"
#include "com_epri_omi_Block.h"

/* Begin Block Signatures
Compiled from Block.java
public class omi.Block extends java.lang.Object 
{
    rowIndexLength     I
    columnIndexLength  I
    valueLength        I
    numRows            I
    numColumns         I
    numElements        I
    offset             Lcom/epri/omi/Base;
    rowIndex           [I
    colIndex           [I
    value              [D
}
End Block Signatures */

/* prototype functions */
void Java_com_epri_omi_block_convertToRowOrder (JNIEnv *, jobject);

#ifndef VERBOSE_INCLUDE
#define VERBOSE_INCLUDE

#endif

/* native methods */

/**
  * Copy scalars from C back to Java and free C version of block.
  * mode = 0,          copy and free
  * mode = OMI_COMMIT, copy w/o free
  * mode = OMI_ABORT,  free w/o copy
  *
  */
void OMI_Block_release
  (JNIEnv *env, jobject jblock, ClassOMI_Block *block, int mode)
{

  jclass         cls;
  jfieldID       fid;
  int VERBOSE = 0;

  /* Check arguments */
  if (env == NULL) {
    printf("Block.c:OMI_Block_release ",
      "!ERROR, NULL 'env' passed as argument.\n");
    return;
  }
  if (jblock == NULL) {
    printf("Block.c:OMI_Block_release ",
      "!ERROR, NULL 'jblock' passed as argument.\n");
    return;
  }
  if (block == NULL) {
    printf("Block.c:OMI_Block_release ",
      "!ERROR, NULL 'block' passed as argument.\n");
    return;
  }
  if ((mode != (jint)NULL) && (mode != OMI_COMMIT) && (mode != OMI_ABORT)) {
    printf("Block.c:OMI_Block_release ",
      "!ERROR, 'mode' is not one of NULL, OMI_COMMIT, or OMI_ABORT.\n");
    return;
  }
  // Announce OMI_Block_release
  if (VERBOSE>0) printf ("OMI_Block_release() Arguments OK.\n");
  if ((mode == (jint)NULL) || (mode == OMI_COMMIT)) {
    // cls
    // In JDK1.1 use the class to get at fields.
    if ((cls = (*env)->GetObjectClass(env, jblock)) == NULL) {
      printf("Block.c:OMI_Block_release ",
        "!ERROR, Cannot GetObjectClass of 'jblock'.\n");
      return;
    }

    // rowIndexLength
    fid = (*env)->GetFieldID(env, cls, "rowIndexLength", "I");
    if (fid == 0) {
      printf ("Block.c:OMI_Block_release ",
        "!ERROR: Cannot GetFieldID of block.rowIndexLength\n");
      return;
    }
    (*env)->SetIntField(env, jblock, fid, block->rowIndexLength); 

    // columnIndexLength
    fid = (*env)->GetFieldID(env, cls, "columnIndexLength", "I");
    if (fid == 0) {
      printf ("Block.c:OMI_Block_release ",
        "!ERROR: Cannot GetFieldID of block.columnIndexLength\n");
      return;
    }
    (*env)->SetIntField(env, jblock, fid, block->columnIndexLength); 

    // valueLength
    fid = (*env)->GetFieldID(env, cls, "valueLength", "I");
    if (fid == 0) {
      printf ("Block.c:OMI_Block_release ",
        "!ERROR: Cannot GetFieldID of block.valueLength\n");
      return;
    }
    (*env)->SetIntField(env, jblock, fid, block->valueLength); 

    // numRows
    fid = (*env)->GetFieldID(env, cls, "numRows", "I");
    if (fid == 0) {
      printf ("Block.c:OMI_Block_release ",
        "!ERROR: Cannot GetFieldID of block.numRows\n");
      return;
    }
    (*env)->SetIntField(env, jblock, fid, block->numRows); 

    // numColumns
    fid = (*env)->GetFieldID(env, cls, "numColumns", "I");
    if (fid == 0) {
      printf ("Block.c:OMI_Block_release ",
        "!ERROR: Cannot GetFieldID of block.numColumns\n");
      return;
    }
    (*env)->SetIntField(env, jblock, fid, block->numColumns); 

    // numElements
    fid = (*env)->GetFieldID(env, cls, "numElements", "I");
    if (fid == 0) {
      printf ("Block.c:OMI_Block_release ",
        "!ERROR: Cannot GetFieldID of block.numElements\n");
      return;
    }
    (*env)->SetIntField(env, jblock, fid, block->numElements); 

    // blockType
    fid = (*env)->GetFieldID(env, cls, "blockType", "I");
    if (fid == 0) {
      printf ("Block.c:OMI_Block_release ",
        "!ERROR: Cannot GetFieldID of block.blockType\n");
      return;
    }
    (*env)->SetIntField(env, jblock, fid, block->blockType); 

  }  // if NULL or OMI_COMMIT
  
  if ((mode == (jint)NULL) ||  (mode == OMI_ABORT)) {
  // copy back elements and free C array
    (*env)->ReleaseIntArrayElements(env, block->jrowIndex, block->rowIndex, 0);
    (*env)->ReleaseIntArrayElements(env, block->jcolIndex, block->colIndex, 0);
    (*env)->ReleaseDoubleArrayElements(env, block->jvalue, block->value, 0);
    free(block);
  }  // if NULL or OMI_ABORT
    
} // void OMI_Block_release()

// ------------------------------------------------------------------
ClassOMI_Block *OMI_Block_get(JNIEnv *env, jobject jblock)
{ 

  jclass          cls;
  jfieldID        fid;
  ClassOMI_Block  *block;
  
  // Check arguments.
  if (env == NULL) {
    printf("Block.c:OMI_Block_get ",
      "!ERROR, NULL 'env' passed as argument.\n");
    return NULL;
  }  
  if (jblock == NULL) {
    printf("Block.c:OMI_Block_get ",
      "!ERROR, NULL 'jblock' passed as argument.\n");
    return NULL;
  }

  // Announce OMI_Block_get
//	printf ("Block->OMI_Block_get() Attributes OK.\n"); 
  // malloc space for a block
  if ((block = (ClassOMI_Block *) malloc(sizeof(ClassOMI_Block))) == NULL) {
    printf("Model.c:OMI_Block_get ",
      "!ERROR, Cannot malloc space for block\n");
    return NULL;
  }
  
  // cls
  // In JDK1.1 use the class to get at fields.
  if ((cls = (*env)->GetObjectClass(env, jblock)) == NULL) {
    printf("Block.c:OMI_Block_get ",
      "!ERROR, Cannot GetObjectClass of 'jblock'.\n");
    OMI_Block_release(env, jblock, block, OMI_ABORT);
    return NULL;
  }

  // load local array references and lengths
  // rowIndex
  fid = (*env)->GetFieldID(env, cls, "rowIndex", "[I");
  if (fid == 0) {
    printf ("Block.c:OMI_Block_get ",
      "!ERROR: Cannot GetFieldID of 'rowIndex'\n");
    return NULL;
  }
  if ((block->jrowIndex = (jintArray) (*env)->GetObjectField(env,jblock,fid)) == NULL) {
    printf("Block.c:OMI_Block_get ",
      "!ERROR, Cannot GetObjectField of 'rowIndex'.\n");
    return NULL;
  }
  if ((block->rowIndex = (*env)->GetIntArrayElements(env, block->jrowIndex, 0)) == NULL) {
    printf("Block.c:OMI_Block_get ",
      "!ERROR, Cannot GetIntArrayElements of 'block->jrowIndex'.\n");
    return NULL;
  }
  block->rowIndexLength = (*env)->GetArrayLength(env, block->jrowIndex);  // Exception? TBD
  
  // colIndex
  fid = (*env)->GetFieldID(env, cls, "colIndex", "[I");
  if (fid == 0) {
    printf ("Block.c:OMI_Block_get ",
      "!ERROR: Cannot GetFieldID of 'colIndex'\n");
    return NULL;
  }
  if ((block->jcolIndex = (jintArray) (*env)->GetObjectField(env,jblock,fid)) == NULL) {
    printf("Block.c:OMI_Block_get ",
      "!ERROR, Cannot GetObjectField of 'colIndex'.\n");
    return NULL;
  }
  if ((block->colIndex = (*env)->GetIntArrayElements(env, block->jcolIndex, 0)) == NULL) {
    printf("Block.c:OMI_Block_get ",
      "!ERROR, Cannot GetIntArrayElements of 'block->jcolIndex'.\n");
    return NULL;
  }
  block->columnIndexLength = (*env)->GetArrayLength(env, block->jcolIndex);
  
  // value
  fid = (*env)->GetFieldID(env, cls, "value", "[D");
  if (fid == 0) {
    printf ("Block.c:OMI_Block_get ",
    "!ERROR: Cannot GetFieldID of 'value'\n");
    // TBD Here we release value w/o a GetDoubleArrayElements
    return NULL;
  }
  if ((block->jvalue = (jintArray) (*env)->GetObjectField(env,jblock,fid)) == NULL) {
    printf("Block.c:OMI_Block_get ",
      "!ERROR, Cannot GetObjectField of 'value'.\n");
    return NULL;
  }
  if ((block->value = (*env)->GetDoubleArrayElements(env, block->jvalue, 0)) == NULL) {
    printf("Block.c:OMI_Block_get ",
      "!ERROR, Cannot GetDoubleArrayElements of 'block->jvalue'.\n");
    return NULL;
  }
  block->valueLength = (*env)->GetArrayLength(env, block->jvalue);
  //
  // numRows is treated special as the number of actual rows in the Block
  fid = (*env)->GetFieldID(env, cls, "numRows", "I");
  if (fid == 0) {
    printf ("Block.c:OMI_Block_get !ERROR: Cannot GetFieldID of 'numRows'\n");
    return NULL;
  }
  block->numRows = (*env)->GetIntField(env,jblock,fid);

  //
  // numColumns is treated special as the number of actual columns in the block
  fid = (*env)->GetFieldID(env, cls, "numColumns", "I");
  if (fid == 0) {
    printf ("Block.c:OMI_Block_get !ERROR: Cannot GetFieldID of 'numColumns'\n");
    return NULL;
  }
  block->numColumns = (*env)->GetIntField(env,jblock,fid);

  //
  // numElements is treated special as the number of actual
  // values in the value field
  fid = (*env)->GetFieldID(env, cls, "numElements", "I");
  if (fid == 0) {
    printf ("Block.c:OMI_Block_get !ERROR: Cannot GetFieldID of 'numElements'\n");
    return NULL;
  }
  block->numElements = (*env)->GetIntField(env,jblock,fid);

  // blockType
  fid = (*env)->GetFieldID(env, cls, "blockType", "I");
  if (fid == 0) {
    printf ("Block.c:OMI_Block_get ",
      "!ERROR: Cannot GetFieldID of 'blockType'\n");
    return NULL;
  }
  block->blockType = (*env)->GetIntField(env,jblock,fid);

  // Get the offset
  fid = (*env)->GetFieldID(env, cls, "offset", "Lcom/epri/omi/Base;");
  if (fid == 0) {
    printf ("Block.c:OMI_Block_get ",
      "!ERROR: Cannot GetFieldID of 'offset'\n");
    return NULL;
  }
  if ((block->joffset = (*env)->GetObjectField(env, jblock, fid)) == NULL) {
    printf ("Block.c:OMI_Block_get ",
      "!ERROR: Cannot GetObjectField of 'offset'\n");
    return NULL;
  }

  // success!
  return block;
} //  ClassOMI_Block OMI_Block_get()

// ------------------------------------------------------------------
void OMI_Block_oprint
  (JNIEnv *env, jobject jblock, jstring jprefix, FILE *out)
{ 
/*
  extern strlen(char *);
  extern strcat(char *, char *);
  extern strcpy(char *, int);
*/ 
  const char      *prefix;
  jstring         jnewPrefix;
  ClassOMI_Block  *block;
  
  integer i;
  char newPrefix[OMI_MaxStrLen], 
       offsetName[OMI_MaxStrLen];

  // Check arguments.
  if (env == NULL) {
    fprintf(out, "Block.c:OMI_Block_oprint ",
      "!ERROR, NULL 'env' passed as argument.\n");
    return;
  }  
  if (jblock == NULL) {
    fprintf(out, "Block.c:OMI_Block_oprint ",
      "!ERROR, NULL 'jblock' passed as argument.\n");
    return;
  }
  if (jprefix == NULL) {
    fprintf(out, "Block.c:OMI_Block_oprint ",
      "!ERROR, NULL 'jprefix' passed as argument.\n");
    return;
  }

  // Announce oprint
  // In JDK 1.1 get a local reference to a string.
  // The UTF string is like a C string.
  // The alternative is a unicode string.
//JDK102 	strcpy(prefix, makeCString(p));
  if ((prefix = (*env)->GetStringUTFChars(env, jprefix, 0)) == NULL) {
    fprintf(out, "Block.c:OMI_Block_oprint ",
      "!ERROR, Cannot GetStringUTFChars of 'jprefix'.\n");
    return;
  }
//	strcpy(prefix, prefix);
  fprintf(out,"%sBlock->oprint() Attributes OK.\n", prefix);
   
  // get block
  if ((block = OMI_Block_get(env, jblock)) == NULL) {
    fprintf(out, "Block.c:OMI_Block_oprint "\
      "!ERROR, Cannot OMI_Block_get of 'block'.\n");
    return;
  }

  // print block dimensions
  fprintf(out, "%sBlock->rowIndex.length @ %X = %d\n",
    prefix, block->rowIndex, block->rowIndexLength);
  fprintf(out, "%sBlock->colIndex.length @ %X = %d\n",
    prefix, block->colIndex, block->columnIndexLength);
  fprintf(out, "%sBlock->value.length @ %X = %d\n",
    prefix, block->value, block->valueLength);
  fprintf(out, "%sBlock->numRows @ %X = %d\n",
    prefix, &block->numRows, block->numRows);
  fprintf(out, "%sBlock->numColumns @ %X = %d\n",
    prefix, &block->numColumns, block->numColumns);
  fprintf(out, "%sBlock->numElements @ %X = %d\n",
    prefix, &block->numElements, block->numElements);
  fprintf(out, "%sBlock->blockType @ %X = %d\n",
    prefix, &block->blockType, block->blockType);

  // copy prefix into newPrefix, then concat Block->offset;
  sprintf(offsetName, "Block->offset @ %X, ", &block->offset);  // Exception? TBD
  if (strlen(offsetName) + strlen(prefix) > OMI_MaxStrLen) {
    fprintf(out, "Block.c:OMI_Block_oprint ",
      "!ERROR: Cannot strcat,  strlen(newPrefix) (= %d) > %d\n",
      strlen(offsetName) + strlen(prefix), OMI_MaxStrLen);
    fprintf(out, "  offsetName = '%s'\n", offsetName       );
    fprintf(out, "  prefix     = '%s'\n", prefix       );
    return;
  }
  strcat(strcpy(newPrefix,prefix),offsetName);  
  if ((jnewPrefix = (jstring) (*env)->NewStringUTF(env, newPrefix)) == NULL) {
    fprintf(out, "Block.c:OMI_Block_oprint ",
      "!ERROR: Cannot NewStringUTF of 'newPrefix'\n");
    return;
  }
  
  OMI_Base_oprint(env, block->joffset, jnewPrefix, out);
  // Do not seem to need to release strings from NewStringUTF
  
  fprintf(out, "%sBlock->rowindex @ %X =\n  ", prefix, block->rowIndex);
  for (i=0; i < block->rowIndexLength; i++) fprintf(out, "%d ", block->rowIndex[i]); 
  fprintf(out, "\n");
  fprintf(out, "%sBlock->colindex @ %X =\n  ", prefix, block->colIndex);
  for (i=0; i < block->columnIndexLength; i++) fprintf(out, "%d ", block->colIndex[i]); 
  fprintf(out, "\n");
  fprintf(out, "%sBlock->value    @ %X =\n  ", prefix, block->value);
  for (i=0; i < block->valueLength; i++) fprintf(out, "%g ", block->value[i]   ); 
  fprintf(out, "\n");

  // release booty gotten
  OMI_Block_release(env, jblock, block, (jint)NULL);
  (*env)->ReleaseStringUTFChars(env, jprefix, prefix);
} // void OMI_Block_oprint()

#pragma export on
// ------------------------------------------------------------------

/**
  * Class:     com_epri_omi_Block
  * Method:    pprint
  * Signature: (Ljava/lang/String;)V
  */
JNIEXPORT void JNICALL Java_com_epri_omi_Block_pprint
  (JNIEnv *env, jobject jblock, jstring jprefix)
{
  FILE *out;
  time_t  t = time(NULL);
  if ((out = fopen("Java_Native_Output", "a")) == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_pprint ",
      "!ERROR, Cannot fopen Java_Native_Output.\n");
    return;
  }
  fprintf(out,"\nBlock->pprint('') %s", asctime(localtime(&t)));
  OMI_Block_oprint(env, jblock, jprefix, out);
  fclose(out);
}  // void JNICALL Java_com_epri_omi_Block_pprint()

// ------------------------------------------------------------------
/*JDK102
integer OMI_Block_get_rowIndex_length(struct HOMI_Block *o)
{ 
  struct ClassOMI_Block  *block    = (struct ClassOMI_Block  *) unhand(o);
  printf("Block.c:Java_com_epri_omi_Block_get_rowIndex_length: returns block->rowIndexLength = %d\n", 
    block->rowIndexLength);
  return block->rowIndexLength;
}
JDK102*/
integer Java_com_epri_omi_Block_get_rowIndex_length(JNIEnv *env, jobject jblock)
{ 
  jclass     cls;
  jfieldID   fid;
  jintArray  jrowIndex;
  
  // Check arguments.
  if (env == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_get_rowIndex_length ",
      "!ERROR, NULL 'env' passed as argument.\n");
    return -1;
  }
  if (jblock == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_get_rowIndex_length ",
      "!ERROR, NULL 'jblock' passed as argument.\n");
    return -1;
  }

//JDK102	struct ClassOMI_Block  *block    = (struct ClassOMI_Block  *) unhand(o);
  // In JDK1.1 use the class to get at fields.
  if ((cls = (*env)->GetObjectClass(env, jblock)) == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_get_rowIndex_length ",
      "!ERROR, Cannot GetObjectClass of 'jblock'.\n");
    return -1;
  }
  fid = (*env)->GetFieldID(env, cls, "rowIndex", "[I");
  if (fid == 0) {
    printf ("Block.c:Java_com_epri_omi_Block_get_rowIndex_length ",
      "!ERROR: Cannot GetFieldID of 'rowIndex'\n");
    return -1;
  }
  if ((jrowIndex = (jintArray) (*env)->GetObjectField(env, jblock, fid)) == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_get_rowIndex_length ",
      "!ERROR, Cannot GetObjectField of 'rowIndex'.\n");
    return -1;
  }
  /*
  printf("Java_com_epri_omi_Block_get_rowIndex_length: returns block->rowIndexLength = %d\n", 
    (*env)->GetArrayLength(env, jrowIndex));
  */
  return (*env)->GetArrayLength(env, jrowIndex);
} //  integer Java_com_epri_omi_Block_get_rowIndex_length()

// ------------------------------------------------------------------
/*JDK102
integer OMI_Block_get_colIndex_length(struct HOMI_Block *o)
{ 
  struct ClassOMI_Block  *block    = (struct ClassOMI_Block  *) unhand(o);
  printf("Block.c:Java_com_epri_omi_Block_get_colIndex_length: returns block->columnIndexLength = %d\n", 
    block->columnIndexLength);
  return block->columnIndexLength;
}
JDK102*/
integer Java_com_epri_omi_Block_get_colIndex_length(JNIEnv *env, jobject jblock)
{ 
  jclass     cls;
  jfieldID   fid;
  jintArray  jcolIndex;
  
  // Check arguments.
  if (env == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_get_colIndex_length ",
      "!ERROR, NULL 'env' passed as argument.\n");
    return -1;
  }
  if (jblock == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_get_colIndex_length ",
      "!ERROR, NULL 'jblock' passed as argument.\n");
    return -1;
  }

  // In JDK1.1 use the class to get at fields.
  if ((cls = (*env)->GetObjectClass(env, jblock)) == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_get_colIndex_length ",
      "!ERROR, Cannot GetObjectClass of 'jblock'.\n");
    return -1;
  }
  fid = (*env)->GetFieldID(env, cls, "colIndex", "[I");
  if (fid == 0) {
    printf ("Block.c:Java_com_epri_omi_Block_get_colIndex_length ",
      "!ERROR: Cannot GetFieldID of 'colIndex'\n");
    return -1;
  }
  if ((jcolIndex = (jintArray) (*env)->GetObjectField(env, jblock, fid)) == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_get_colIndex_length ",
      "!ERROR, Cannot GetObjectField of 'colIndex'.\n");
    return -1;
  }
  /*
  printf("Java_com_epri_omi_Block_get_colIndex_length: returns block->columnIndexLength = %d\n", 
    (*env)->GetArrayLength(env, jcolIndex));
  */
  return (*env)->GetArrayLength(env, jcolIndex);
} //  integer Java_com_epri_omi_Block_get_colIndex_length()

// ------------------------------------------------------------------
/*JDK102
integer OMI_Block_get_value_length(struct HOMI_Block *o)
{ 
  struct ClassOMI_Block  *block    = (struct ClassOMI_Block  *) unhand(o);
  printf("Block.c:Java_com_epri_omi_Block_get_value_length: returns block->valueLength = %d\n", 
    block->valueLength);
  return block->valueLength;
}
JDK102*/
integer Java_com_epri_omi_Block_get_value_length(JNIEnv *env, jobject jblock)
{ 
  jclass        cls;
  jfieldID      fid;
  jdoubleArray  jvalue;
  
  // Check arguments.
  if (env == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_get_value_length ",
      "!ERROR, NULL 'env' passed as argument.\n");
    return -1;
  }
  if (jblock == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_get_value_length ",
      "!ERROR, NULL 'jblock' passed as argument.\n");
    return -1;
  }

  // In JDK1.1 use the class to get at fields.
  if ((cls = (*env)->GetObjectClass(env, jblock)) == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_get_value_length ",
      "!ERROR, Cannot GetObjectClass of 'jblock'.\n");
    return -1;
  }
  fid = (*env)->GetFieldID(env, cls, "value", "[D");
  if (fid == 0) {
    printf ("Block.c:Java_com_epri_omi_Block_get_value_length ",
      "!ERROR: Cannot GetFieldID of 'value'\n");
    return -1;
  }
  if ((jvalue = (jdoubleArray) (*env)->GetObjectField(env,jblock,fid)) == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_get_value_length ",
      "!ERROR, Cannot GetObjectField of 'value'.\n");
    return -1;
  }
  /*
  printf("Java_com_epri_omi_Block_get_value_length: returns block->valueLength = %d\n", 
    (*env)->GetArrayLength(env, jvalue));
  */
  return (*env)->GetArrayLength(env, jvalue);
} //  integer Java_com_epri_omi_Block_get_value_length()

// ------------------------------------------------------------------
/*JDK102
integer OMI_Block_get_num_elems(struct HOMI_Block *o)
{ 
  struct ClassOMI_Block  *block    = (struct ClassOMI_Block  *) unhand(o);
  printf("Block.c:Java_com_epri_omi_Block_get_rowIndex: returns block->numElements = %d\n", block->numElements);
  return block->numElements;
}
JDK102*/
integer Java_com_epri_omi_Block_get_num_elems(JNIEnv *env, jobject jblock)
{ 
  jclass     cls;
  jfieldID   fid;
  jint       jnum_elems;
  
  // Check arguments.
  if (env == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_get_num_elems ",
      "!ERROR, NULL 'env' passed as argument.\n");
    return -1;
  }  
  if (jblock == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_get_num_elems ",
      "!ERROR, NULL 'jblock' passed as argument.\n");
    return -1;
  }

  // In JDK1.1 use the class to get at fields.
  if ((cls = (*env)->GetObjectClass(env, jblock)) == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_get_num_elems ",
      "!ERROR, Cannot GetObjectClass of 'jblock'.\n");
    return -1;
  }
  fid = (*env)->GetFieldID(env, cls, "numElements", "I");
  if (fid == 0) {
    printf ("Block.c:Java_com_epri_omi_Block_get_num_elems ",
      "!ERROR: Cannot GetFieldID of 'value'\n");
    return -1;
  }
  jnum_elems = (jint) (*env)->GetIntField(env,jblock,fid);  // Exception? TBD
  /*
  printf("Java_com_epri_omi_Block_get_num_elems: returns block->numElements = %d\n", 
    jnum_elems);
  */
  return jnum_elems;
} //  integer Java_com_epri_omi_Block_get_num_elems()

// ------------------------------------------------------------------
/*JDK102
integer *OMI_Block_get_rowIndex (struct HOMI_Block *o)
{ 
  struct ClassOMI_Block  *block    = (struct ClassOMI_Block  *) unhand(o);
  integer *result = (integer *) unhand(block->rowIndex);
  printf("Block.c:Java_com_epri_omi_Block_get_rowIndex: block->rowIndex @ %X, *result = %d\n",
    result, *result);
  return (integer *) unhand(block->rowIndex);
}
JDK102*/

// ------------------------------------------------------------------
integer *Java_com_epri_omi_Block_get_rowIndex(JNIEnv *env, jobject jblock)
{ 
  jclass     cls;
  jfieldID   fid;
  jintArray  jrowIndex;
  jint       *rowIndex;
  
  // Check arguments.
  if (env == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_get_rowIndex ",
      "!ERROR, NULL 'env' passed as argument.\n");
    return NULL;
  }  
  if (jblock == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_get_rowIndex ",
      "!ERROR, NULL 'jblock' passed as argument.\n");
    return NULL;
  }

  // In JDK1.1 use the class to get at fields.
  if ((cls = (*env)->GetObjectClass(env, jblock)) == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_get_rowIndex ",
      "!ERROR, Cannot GetObjectClass of 'jblock'.\n");
    return NULL;
  }
  fid = (*env)->GetFieldID(env, cls, "rowIndex", "[I");
  if (fid == 0) {
    printf ("Block.c:Java_com_epri_omi_Block_get_rowIndex ",
      "!ERROR: Cannot GetFieldID of 'rowIndex'\n");
    return NULL;
  }
  if ((jrowIndex = (jintArray) (*env)->GetObjectField(env, jblock, fid)) == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_get_rowIndex ",
      "!ERROR, Cannot GetObjectField of 'rowIndex'.\n");
    return NULL;
  }
  if ((rowIndex = (*env)->GetIntArrayElements(env, jrowIndex, 0)) == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_get_rowIndex ",
      "!ERROR, Cannot GetIntArrayElements of 'jrowIndex'.\n");
    return NULL;
  }
  /*
  printf("Block.c:Java_com_epri_omi_Block_get_rowIndex: block->rowIndex @ %X, *rowIndex = %d\n",
    rowIndex, *rowIndex);
  */
  return (integer *)rowIndex;
} //integer *Java_com_epri_omi_Block_get_rowIndex()

// ------------------------------------------------------------------
/*JDK102
integer *OMI_Block_get_colIndex (struct HOMI_Block *o)
{ 
  struct ClassOMI_Block  *block    = (struct ClassOMI_Block  *) unhand(o);
  integer *result = (integer *) unhand(block->colIndex);
  printf("Block.c:Java_com_epri_omi_Block_get_colIndex: block->colIndex @ %X, *result = %d\n",
    result, *result);
  return (integer *) unhand(block->colIndex);
}
JDK102*/
integer *Java_com_epri_omi_Block_get_colIndex(JNIEnv *env, jobject jblock)
{ 
  jclass     cls;
  jfieldID   fid;
  jintArray  jcolIndex;
  jint       *colIndex;
  
  // Check arguments.
  if (env == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_get_colIndex ",
      "!ERROR, NULL 'env' passed as argument.\n");
    return NULL;
  }  
  if (jblock == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_get_colIndex ",
      "!ERROR, NULL 'jblock' passed as argument.\n");
    return NULL;
  }

  // In JDK1.1 use the class to get at fields.
  if ((cls = (*env)->GetObjectClass(env, jblock)) == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_get_colIndex ",
      "!ERROR, Cannot GetObjectClass of 'jblock'.\n");
    return NULL;
  }
  fid = (*env)->GetFieldID(env, cls, "colIndex", "[I");
  if (fid == 0) {
    printf ("Block.c:Java_com_epri_omi_Block_get_colIndex ",
      "!ERROR: Cannot GetFieldID of 'colIndex'\n");
    return NULL;
  }
  if ((jcolIndex = (jintArray) (*env)->GetObjectField(env, jblock, fid)) == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_get_colIndex ",
      "!ERROR, Cannot GetObjectField of 'colIndex'.\n");
    return NULL;
  }
  if ((colIndex = (*env)->GetIntArrayElements(env, jcolIndex, 0)) == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_get_colIndex ",
      "!ERROR, Cannot GetIntArrayElements of 'jcolIndex'.\n");
    return NULL;
  }
  /*
  printf("Block.c:Java_com_epri_omi_Block_get_colIndex: block->colIndex @ %X, *colIndex = %d\n",
    colIndex, *colIndex);
  */
  return (integer *)colIndex;
} //  integer *Java_com_epri_omi_Block_get_colIndex()

// ------------------------------------------------------------------
/*JDK102
doublereal *OMI_Block_get_value(struct HOMI_Block *o)
{ 
  struct ClassOMI_Block  *block    = (struct ClassOMI_Block  *) unhand(o);
  doublereal *result = (doublereal *) unhand(block->value);
  printf("Block.c:Java_com_epri_omi_Block_get_value: block->value @ %X, *result = %g\n",
    result, *result);
  return (doublereal *) unhand(block->value);
}
JDK102*/
doublereal *Java_com_epri_omi_Block_get_value(JNIEnv *env, jobject jblock)
{ 
  jclass        cls;
  jfieldID      fid;
  jdoubleArray  jvalue;
  jdouble       *value;
  
  // Check arguments.
  if (env == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_get_value ",
      "!ERROR, NULL 'env' passed as argument.\n");
    return NULL;
  }  
  if (jblock == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_get_value ",
      "!ERROR, NULL 'jblock' passed as argument.\n");
    return NULL;
  }

  // In JDK1.1 use the class to get at fields.
  if ((cls = (*env)->GetObjectClass(env, jblock)) == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_get_value ",
      "!ERROR, Cannot GetObjectClass of 'jblock'.\n");
    return NULL;
  }
  fid = (*env)->GetFieldID(env, cls, "value", "[D");
  if (fid == 0) {
    printf ("Block.c:Java_com_epri_omi_Block_get_value ",
      "!ERROR: Cannot GetFieldID of 'value'\n");
    return NULL;
  }
  if ((jvalue = (jdoubleArray) (*env)->GetObjectField(env, jblock, fid)) == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_get_value ",
      "!ERROR, Cannot GetObjectField of 'value'.\n");
    return NULL;
  }
  if ((value = (*env)->GetDoubleArrayElements(env, jvalue, 0)) == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_get_value ",
      "!ERROR, Cannot GetIntArrayElements of 'jvalue'.\n");
    return NULL;
  }
  /*
  printf("Block.c:Java_com_epri_omi_Block_get_value: block->value @ %X, *value = %d\n",
    value, *value);
  */
  return value;
} //  doublereal *Java_com_epri_omi_Block_get_value()

// ------------------------------------------------------------------
/*JDK102
void OMI_Block_set_num_elems(struct HOMI_Block *o, integer len)
{ 
  struct ClassOMI_Block  *block    = (struct ClassOMI_Block  *) unhand(o);
  block->numElements = len;
  printf("Block.c:Java_com_epri_omi_Block_set_num_elems: block->numElements @ %X = %d\n",
    &block->numElements, block->numElements);
}
JDK102*/
integer Java_com_epri_omi_Block_set_num_elems(JNIEnv *env, jobject jblock, jint len)
{ 
  jclass     cls;
  jfieldID   fid;
  jint       jnum_elems;
  FILE       *out;
  time_t     t;
  
  // Check arguments.
  if (env == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_set_num_elems ",
      "!ERROR, NULL 'env' passed as argument.\n");
    return 0;
  }  
  if (jblock == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_set_num_elems ",
      "!ERROR, NULL 'jblock' passed as argument.\n");
    return 0;
  }
  if (len < 0) {
    printf("Block.c:Java_com_epri_omi_Block_set_num_elems ",
      "!ERROR, Arguement 'len' (= %d) < 0.\n");
    return 0;
  }

  // In JDK1.1 use the class to get at fields.
  if ((cls = (*env)->GetObjectClass(env, jblock)) == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_set_num_elems ",
      "!ERROR, Cannot GetObjectClass of 'jblock'.\n");
    return 0;
  }
  fid = (*env)->GetFieldID(env, cls, "numElements", "I");
  if (fid == 0) {
    printf ("Block.c:Java_com_epri_omi_Block_set_num_elems ",
      "!ERROR: Cannot GetFieldID of 'value'\n");
    return 0;
  }
  (*env)->SetIntField(env, jblock, fid, len);  // Exception? TBD
  
  if (1) {
    // Open output file for append.
    if ((out = fopen("Java_Native_Output", "a")) == NULL) {
      fprintf(out,"Block.c:Java_com_epri_omi_Block_set_num_elems ",
        "!ERROR, Cannot fopen Java_Native_Output.\n");
      return (integer)NULL;
    }
    
    // print entry time
    t = time(NULL);
    fprintf(out, "\nBlock.c:Java_com_epri_omi_Block_set_num_elems(filename) %s", asctime(localtime(&t)));

    // print set value
    jnum_elems = (*env)->GetIntField(env, jblock, fid);
    fprintf(out, "Block.c:Java_com_epri_omi_Block_set_num_elems: sets block->numElements = %d\n", 
      jnum_elems);
      
    // close output file
    fprintf(out, "Block.c:Java_com_epri_omi_Block_set_num_elems: return, closing output file\n");
    fclose(out);
  }
  return 1;
} //  integer Java_com_epri_omi_Block_set_num_elems()

//
// ------------------------------------------------------------------
integer Java_com_epri_omi_Block_release_rowIndex
  (JNIEnv *env, jobject jblock, jint *rowIndex)  
{
  jclass     cls;
  jfieldID   fid;
  jintArray  jrowIndex;
  
  // Check arguments.
  if (env == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_release_rowIndex ",
      "!ERROR, NULL 'env' passed as argument.\n");
    return 0;
  }  
  if (jblock == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_release_rowIndex ",
      "!ERROR, NULL 'jblock' passed as argument.\n");
    return 0;
  }
  if (rowIndex == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_release_rowIndex ",
      "!ERROR, NULL 'rowIndex' passed as argument.\n");
    return 0;
  }

  // In JDK1.1 use the class to get at fields.
  if ((cls = (*env)->GetObjectClass(env, jblock)) == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_release_rowIndex ",
      "!ERROR, Cannot GetObjectClass of 'jblock'.\n");
    return 0;
  }
  fid = (*env)->GetFieldID(env, cls, "rowIndex", "[I");
  if (fid == 0) {
    printf ("Block.c:Java_com_epri_omi_Block_release_rowIndex ",
      "!ERROR: Cannot GetFieldID of 'rowIndex'\n");
    return 0;
  }
  if ((jrowIndex = (jintArray) (*env)->GetObjectField(env, jblock, fid)) == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_release_rowIndex ",
      "!ERROR, Cannot GetObjectField of 'rowIndex'.\n");
    return 0;
  }
  // mode = 0, copy back elements and free C array
  (*env)->ReleaseIntArrayElements(env, jrowIndex, rowIndex, 0);
  printf("Java_com_epri_omi_Block_release_rowIndex: block->rowIndex @ %X, *rowIndex = %d\n",
    rowIndex, *rowIndex);
  return 1;
} //  integer Java_com_epri_omi_Block_release_rowIndex

//
// ------------------------------------------------------------------
integer Java_com_epri_omi_Block_release_colIndex
  (JNIEnv *env, jobject jblock, jint *colIndex)  
{
  jclass     cls;
  jfieldID   fid;
  jintArray  jcolIndex;
  
  // Check arguments.
  if (env == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_release_colIndex ",
      "!ERROR, NULL 'env' passed as argument.\n");
    return 0;
  }  
  if (jblock == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_release_colIndex ",
      "!ERROR, NULL 'jblock' passed as argument.\n");
    return 0;
  }
  if (colIndex == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_release_colIndex ",
      "!ERROR, NULL 'colIndex' passed as argument.\n");
    return 0;
  }

  // In JDK1.1 use the class to get at fields.
  if ((cls = (*env)->GetObjectClass(env, jblock)) == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_release_colIndex ",
      "!ERROR, Cannot GetObjectClass of 'jblock'.\n");
    return 0;
  }
  fid = (*env)->GetFieldID(env, cls, "colIndex", "[I");
  if (fid == 0) {
    printf ("Block.c:Java_com_epri_omi_Block_release_colIndex ",
      "!ERROR: Cannot GetFieldID of 'colIndex'\n");
    return 0;
  }
  if ((jcolIndex = (jintArray) (*env)->GetObjectField(env, jblock, fid)) == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_release_colIndex ",
      "!ERROR, Cannot GetObjectField of 'colIndex'.\n");
    return 0;
  }
  // mode = 0, copy back elements and free C array
  (*env)->ReleaseIntArrayElements(env, jcolIndex, colIndex, 0);
  printf("Java_com_epri_omi_Block_release_colIndex: block->colIndex @ %X, *colIndex = %d\n",
    colIndex, *colIndex);
  return 1;
} //  integer Java_com_epri_omi_Block_release_colIndex()

//
// ------------------------------------------------------------------
integer Java_com_epri_omi_Block_release_value
  (JNIEnv *env, jobject jblock, jdouble *value)  
{
  jclass     cls;
  jfieldID   fid;
  jintArray  jvalue;
  
  // Check arguments.
  if (env == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_release_value ",
      "!ERROR, NULL 'env' passed as argument.\n");
    return 0;
  }  
  if (jblock == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_release_value ",
      "!ERROR, NULL 'jblock' passed as argument.\n");
    return 0;
  }
  if (value == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_release_value ",
      "!ERROR, NULL 'value' passed as argument.\n");
    return 0;
  }

  // In JDK1.1 use the class to get at fields.
  if ((cls = (*env)->GetObjectClass(env, jblock)) == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_release_value ",
      "!ERROR, Cannot GetObjectClass of 'jblock'.\n");
    return 0;
  }
  fid = (*env)->GetFieldID(env, cls, "value", "[I");
  if (fid == 0) {
    printf ("Block.c:Java_com_epri_omi_Block_release_value ",
      "!ERROR: Cannot GetFieldID of 'value'\n");
    return 0;
  }
  if ((jvalue = (jdoubleArray) (*env)->GetObjectField(env, jblock, fid)) == NULL) {
    printf("Block.c:Java_com_epri_omi_Block_release_value ",
      "!ERROR, Cannot GetObjectField of 'value'.\n");
    return 0;
  }
  (*env)->ReleaseDoubleArrayElements(env, jvalue, value, 0);
  printf("Java_com_epri_omi_Block_release_value: block->value @ %X, *value = %d\n",
    value, *value);
  return 1;
} //  integer Java_com_epri_omi_Block_release_value()

#pragma export off

void Java_com_epri_omi_block_convertToRowOrder(JNIEnv *env, jobject jblock) {
  const char prefix[] = "Block.c:Java_com_epri_omi_block_convertToRowOrder():";
	const int VERBOSE = 0;
  
  // local variables
  jclass    cls;  // holds pointer to omi.Model.class
  jmethodID mid;  // holds offset for omi.Model.getObjective()

  char javaName[] = "convertToRowOrder";
  char javaSig[]  = "()V"; // for the omi version

  // Check arguments.
  if (env == NULL) {
    printf("%s!ERROR, NULL 'env' passed as argument.\n", prefix);
    //fflush(stdout);
    return;
  };  
  if (jblock == NULL) {
    printf("%s!ERROR, NULL 'jblock' passed as argument.\n", prefix);
    //fflush(stdout);
    return;
  };

  // Announce OMI_Model_release
  if (VERBOSE>0) printf("%sArguments OK.\n",prefix);

  // get the class pointer and method offset
  cls = (*env)->GetObjectClass(env, jblock);
  mid = (*env)->GetMethodID(env, cls, javaName, javaSig);
  if (mid == 0) {
    printf("%s ERROR: Cannot GetMethodID for '%s', '%s'\n",prefix,javaName,javaSig);
    printf("%s return\n",prefix);
    return;
  }
  
  if (VERBOSE>0) // call back to Java
  printf("%s In C, calling Block.convertToRowOrder\n", prefix);
  (*env)->CallVoidMethod(env, jblock, mid);
  if (VERBOSE>0) printf("%s In C, finished Block.convertToRowOrder\n", prefix);
  
  if (VERBOSE>0) printf("%s return\n",prefix);
  return;
} //  void Java_com_epri_omi_block_convertToRowOrder()

void Java_com_epri_omi_block_convertToColumnOrder(JNIEnv *env, jobject jblock) {
  const char prefix[] = "Block.c:Java_com_epri_omi_block_convertToColumnOrder():";
	const int VERBOSE = 0;
  
  // local variables
  jclass    cls;  // holds pointer to omi.Model.class
  jmethodID mid;  // holds offset for omi.Model.getObjective()

  char javaName[] = "convertToColumnOrder";
  char javaSig[]  = "()V"; // for the omi version

  // Check arguments.
  if (env == NULL) {
    printf("%s!ERROR, NULL 'env' passed as argument.\n", prefix);
    //fflush(stdout);
    return;
  };  
  if (jblock == NULL) {
    printf("%s!ERROR, NULL 'jblock' passed as argument.\n", prefix);
    //fflush(stdout);
    return;
  };

  // Announce Java_com_epri_omi_block_convertToColumnOrder
  if (VERBOSE>0) printf("%sArguments OK.\n",prefix);

  // get the class pointer and method offset
  cls = (*env)->GetObjectClass(env, jblock);
  mid = (*env)->GetMethodID(env, cls, javaName, javaSig);
  if (mid == 0) {
    printf("%s ERROR: Cannot GetMethodID for '%s', '%s'\n",prefix,javaName,javaSig);
    printf("%s return\n",prefix);
    return;
  }
  
  if (VERBOSE>0) // call back to Java
  printf("%s In C, calling Block.convertToColumnOrder\n", prefix);
  (*env)->CallVoidMethod(env, jblock, mid);
  if (VERBOSE>0) printf("%s In C, finished Block.convertToColumnOrder\n", prefix);
  
  if (VERBOSE>0) printf("%s return\n",prefix);
  return;
} //  void Java_com_epri_omi_block_convertToColumnOrder()

void Java_com_epri_omi_block_convertToTriples(JNIEnv *env, jobject jblock) {
  int VERBOSE = 0;
	
	// local variables
	jclass    cls;  // holds pointer to omi.Model.class
	jmethodID mid;  // holds offset for omi.Model.getObjective()

  char prefix[] = "Block.c:omi_block_convertToTriples__():";
	char javaName[] = "convertToTriples";
	char javaSig[]  = "()V"; // for the omi version

  // get the class pointer and method offset
  cls = (*env)->GetObjectClass(env, jblock);
  mid = (*env)->GetMethodID(env, cls, javaName, javaSig);
  if (mid == 0) {
    printf("%s ERROR: Cannot GetMethodID for '%s', '%s'\n",prefix,javaName,javaSig);
    printf("%s return\n",prefix);
    return;
  }
  
  if (VERBOSE>0) // call back to Java
  printf("%s In C, calling omi_block_convertToTriples__\n", prefix);
  (*env)->CallVoidMethod(env, jblock, mid);
  if (VERBOSE>0) printf("%s In C, back from omi_block_convertToTriples__\n", prefix);
}  // void Java_com_epri_omi_block_convertToTriples()

