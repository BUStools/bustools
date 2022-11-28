#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <vector>
#include <array>
#include <stdint.h>
#include <zlib.h>

#include "Common.hpp"
#include "BUSData.h"
#include "bustools_compress.h"

size_t pfd_blocksize = 512;

/**
 * @brief Encode `num` using fibonacci encoding into buf, starting at bitpos.
 * @pre the next 128 bits in `buf` are 0 starting from `bitpos`, and wrapped around 192.
 * 		0 <= bit_pos < 3*64 == 192
 *
 * @pre num > 0
 * @post buf now contains the fibo encoding of num from [pre(bitpos); post(bitpos)].
 *		bit_pos has changed, 0 <= bitpos < 3*64 == 192.
 * 		0 or more elements in buf have been concentrated with fibo-encoded numbers and thus written to obuf.
 *
 * @note The largest fibonacci number that fits in a uint64_t is the 91st (0-indexed) fibo number, i.e. the first 92 fibonacci numbers fit in a 64-bit uint.
 * 	Fibonacci encoding appends a single bit as a stop code. Hence, the longest fibonacci encoding uses 93 bits.
 *
 * @param num the number to encode, num > 0
 * @param bufsize the size of the output buffer.
 * @param buf array of `bufsize` uint64_t. num is fibonacci-encoded in buf.
 * @param bitpos the bit position in buf where the next fibo encoding of num starts.
 * @param obuf the ostream buffer to write to.
 * @return bool true iff encoding the current number does not write outside of buf
 */

const auto fibo_begin = fibo64.begin();
const auto fibo_end = fibo64.end();

template <typename BUF_t>
bool fiboEncode(const uint64_t num, const size_t bufsize, BUF_t *buf, size_t &bitpos)
{
	constexpr uint32_t word_size = sizeof(BUF_t) * 8;

	const size_t max_bitpos = bufsize * word_size;
	constexpr BUF_t ONE{1};

	const uint32_t curr_byte_pos = bitpos / word_size;

	uint64_t remainder = num;

	// the ith fibonacci number is the largest fibo not greater than remainder
	auto i = std::upper_bound(fibo_begin, fibo_end, remainder) - 1;

	const uint32_t n_bits = (i - fibo_begin) + 2;

	// Encoding the current number would write out of bounds of buf.
	if (bitpos + n_bits > max_bitpos)
		return false;

	uint32_t next_bit_pos = bitpos + n_bits - 1,
			 bit_offset = next_bit_pos % word_size,
			 buf_offset = (next_bit_pos / word_size) % bufsize;

	// Set the stop bit.
	buf[buf_offset] |= ONE << (word_size - 1 - bit_offset);

	++i;
	while (remainder > 0)
	{
		i = std::upper_bound(fibo_begin, i, remainder) - 1;
		next_bit_pos = bitpos + (i - fibo_begin);
		buf_offset = (next_bit_pos / word_size) % bufsize;
		bit_offset = next_bit_pos % word_size;

		buf[buf_offset] |= ONE << (word_size - 1 - bit_offset);
		remainder -= *i;
	}
	bitpos += n_bits;
	return true;
}

/**
 * @brief pack elem into buf starting at bitpos, using exactly b_bits bits.
 * @pre num is representable using `b_bits` bits.
 * @pre buf has at least one element
 * @pre buf has at least two elements if bitpos + b_bits > 8*sizeof(buf).
 *
 * @param b_bits The number of bits to represent elem.
 * @param elem The number to pack, must be representable with at most `b_bits` bits.
 * @param buf The int64_t array where elem should be packed into.
 * @param bitpos The starting point in bits of where to back elem.
 * @return bool: true iff packing of elem saturates buf[0].
 */
template <typename T>
bool pack_int(
	const uint32_t b_bits,
	T elem,
	T *buf,
	uint32_t &bitpos)
{
	constexpr int32_t dest_wordsize = sizeof(T) * 8;

	// |               |               |
	//    ^             ^
	//  bitpos         dest_wordsize
	int32_t shift = dest_wordsize - bitpos - b_bits;
	T carryover = 0;
	constexpr T ONE{1};

	if (shift < 0)
	{
		// b_bits > (dest_wordsize - bitpos)
		// -> number covers two elements

		carryover = elem & (ONE << -shift) - 1;
		*(buf + 1) = carryover << dest_wordsize + shift;
		elem >>= -shift;
	}

	// shift by max(0, shift)
	*buf |= elem << (shift > 0) * shift;

	bitpos = (dest_wordsize - shift) % dest_wordsize;
	return (shift <= 0);
}

/**
 * @brief Encode a block of integers using NewPFD.
 * @pre BUF has a size of `pfd_block`.size() / 64 * `b_bits` elements.
 *
 * @param pfd_block The numbers to encode.
 * @param index_gaps Output: delta encoded indices of exceptions in `pfd_block`.
 * @param exceptions Output: The most significant bits of each exception.
 * @param b_bits The number of bits to use per int for the block.
 * @param min_element The smallest element in `pfd_block`.
 * @param BUF The buffer to pack the elements of pfd_block into.
 */
template <typename SRC_T, typename DEST_T>
void encode_pfd_block(
	std::vector<SRC_T> &pfd_block,
	std::vector<uint32_t> &index_gaps,
	std::vector<SRC_T> &exceptions,
	const uint32_t b_bits,
	const SRC_T min_element,
	DEST_T *BUF)
{
	index_gaps.clear();
	exceptions.clear();

	SRC_T max_elem_bit_mask = (1ULL << b_bits) - 1;
	uint32_t bitpos = 0,
			 idx = 0,
			 last_ex_idx = 0;

	bool do_increment_pointer = 0;

	// Store the elements in the primary block in pfd_buf using `b_bits` bits each.
	for (SRC_T &elem : pfd_block)
	{
		elem -= min_element;
		if (elem > max_elem_bit_mask)
		{
			// store the overflowing, most significant bits as exceptions.
			exceptions.push_back(elem >> b_bits);
			index_gaps.push_back(idx - last_ex_idx);
			last_ex_idx = idx;

			// store the least b significant bits in the frame.
			elem &= max_elem_bit_mask;
		}

		do_increment_pointer = pack_int<DEST_T>(b_bits, elem, BUF, bitpos);
		BUF += do_increment_pointer;
		++idx;
	}
}

/**
 * @brief Encode a block of `block_size` elements in `pfd_block` and write to of.
 *
 * @param block_size The number of elements in the block.
 * @param pfd_block The elements to encode.
 * @param index_gaps Vector to store delta-encoded indices of exceptions.
 * @param exceptions Vector to store the most significant bits of exceptions.
 * @param fibonacci_buf array used for storing fibonacci encoding.
 * @param b_bits The number of bits each num in `pfd_block` is represented with.
 * @param min_element The smallest element in `pfd_block`.
 * @param of The ostream to write out the encoding.
 */

template <typename SRC_T, typename DEST_T>
size_t new_pfd(
	const size_t block_size,
	std::vector<SRC_T> &pfd_block,
	std::vector<uint32_t> &index_gaps,
	std::vector<SRC_T> &exceptions,
	DEST_T *fibonacci_buf,
	const size_t b_bits,
	const SRC_T min_element,
	size_t fibonacci_bufsize,
	DEST_T *PFD_buf)
{
	bool success = true;
	constexpr size_t wordsize = sizeof(DEST_T) * 8;

	size_t buf_size = (block_size * b_bits) / wordsize;
	std::fill(PFD_buf, PFD_buf + buf_size, 0ULL);

	size_t bitpos{0};

	encode_pfd_block<SRC_T, DEST_T>(pfd_block, index_gaps, exceptions, b_bits, min_element, PFD_buf);

	size_t n_exceptions = index_gaps.size();
	// For more compact fibonacci encoding, we pack the fibo encoded values together
	success &= fiboEncode(b_bits + 1, fibonacci_bufsize, fibonacci_buf, bitpos);
	success &= fiboEncode(min_element + 1, fibonacci_bufsize, fibonacci_buf, bitpos);
	success &= fiboEncode(n_exceptions + 1, fibonacci_bufsize, fibonacci_buf, bitpos);

	for (const auto &i : index_gaps)
	{
		// we must increment by one, since the first index gap can be zero
		success &= fiboEncode(i + 1, fibonacci_bufsize, fibonacci_buf, bitpos);
	}
	for (const auto &ex : exceptions)
	{
		// These are always > 0 since they contain the most significant bits of exception
		success &= fiboEncode(ex, fibonacci_bufsize, fibonacci_buf, bitpos);
	}

	size_t n_elems = bitpos / wordsize + (bitpos % wordsize > 0);
	success &= (n_elems + buf_size <= fibonacci_bufsize);

	std::memcpy(fibonacci_buf + n_elems, PFD_buf, buf_size * sizeof(DEST_T) * success);
	n_elems += buf_size;

	return n_elems * success;
}

/**
 * @brief Compute the smallest element, and the number of bits required to encode at least 90% of
 * the numbers in pfd_scratch.
 * @pre pfd_scratch contains the numbers to encode in a block
 * @post pfd_scratch may have been reordered.
 * @param block_size The number of elements in pfd_scratch, i.e. the block size.
 * @param pfd_scratch The numbers to encode in a single block using NewPFD.
 * @param min_element Output: The smallest element in `pfd_scratch`.
 * @param b_bits Output: The minimum number of bits that are enough to encode at least ~90% of the elements in `pfd_scratch`.
 */

template<typename T>
void compute_pfd_params(
	const size_t block_size,
	std::vector<T> &pfd_scratch,
	const T &min_element,
	uint32_t &b_bits)
{
	const size_t nth_elem_idx = (size_t)((double)block_size * 0.9);

	std::nth_element(pfd_scratch.begin(), pfd_scratch.begin() + nth_elem_idx, pfd_scratch.end());
	T nth_max_element = pfd_scratch[nth_elem_idx];

	b_bits = sizeof(T) * 8 - 1 - __builtin_clrsb(nth_max_element - min_element);
	b_bits = b_bits ?: 1;
}

/**
 * @brief Compress barcodes of rows using delta-runlen(0)-fibonacci encoding and write to `of`.
 *
 * @param rows BUSData array, contains at least `row_count` elements
 * @param row_count The number of barcodes to compress.
 * @param of The ostream for writing the encoding to.
 * @return bool true iff encoding does not go out of bounds of obuf
 */
template <typename T>
bool compress_barcodes(BUSData const *const rows, const int row_count, char *obuf, const size_t &obuf_size, size_t &global_bufpos)
{
	bool success = true;
	uint64_t barcode = 0,
			 last_bc = 0,
			 runlen = 0;
	size_t wordsize = sizeof(T) * 8;

	const size_t fibonacci_bufsize = (obuf_size - global_bufpos) / sizeof(T);
	T *fibonacci_buf = (T *)(obuf + global_bufpos);
	size_t bitpos{0};

	for (int i = 0; i < row_count && success; ++i)
	{
		barcode = rows[i].barcode;

		// delta encoding
		barcode -= last_bc;

		// Runlength encoding of zeros
		if (barcode == 0)
		{
			++runlen;
		}
		else
		{
			// Increment values as fibo cannot encode 0
			if (runlen)
			{
				success &= fiboEncode(1ULL, fibonacci_bufsize, fibonacci_buf, bitpos);
				success &= fiboEncode(runlen, fibonacci_bufsize, fibonacci_buf, bitpos);
				runlen = 0;
			}

			if (barcode + 1 == 0)
			{
				std::cerr << "This input file needs sorting. Please sort this file and try again." << std::endl;
				throw std::runtime_error("Input needs sorting prior to compression");
			}

			success &= fiboEncode(barcode + 1, fibonacci_bufsize, fibonacci_buf, bitpos);
		}
		last_bc = rows[i].barcode;
	}

	// Take care of the last run of zeros in the delta-encoded barcodes
	if (runlen)
	{
		success &= fiboEncode(1ULL, fibonacci_bufsize, fibonacci_buf, bitpos);
		success &= fiboEncode(runlen, fibonacci_bufsize, fibonacci_buf, bitpos);
	}

	global_bufpos += (bitpos / wordsize + (bitpos % wordsize > 0)) * sizeof(T);
	return success;
}

/**
 * @brief Compress UMIs of rows using periodic_delta-runlen(0)-fibonacci encoding and write to `of`.
 *
 * @param rows BUSData array, contains at least `row_count` elements
 * @param row_count The number of UMIs to compress.
 * @param of The ostream for writing the encoding to.
 * @return bool true iff encoding does not go out of bounds of obuf
 */
template <typename T>
bool lossless_compress_umis(BUSData const *const rows, const int row_count, char *obuf, const size_t &obuf_size, size_t &global_bufpos)
{
	bool success = true;
	uint64_t last_bc = rows[0].barcode + 1,
			 last_umi = 0,
			 bc, umi, diff;

	size_t wordsize = sizeof(T) * 8;

	const size_t fibonacci_bufsize = (obuf_size - global_bufpos) / sizeof(T);
	T *fibonacci_buf = (T *)(obuf + global_bufpos);
	size_t bitpos{0};

	uint64_t runlen{0};
	const uint32_t RLE_val{0ULL};

	for (int i = 0; i < row_count && success; ++i)
	{
		bc = rows[i].barcode;

		// We must increment umi, since a UMI==0 will confuse the runlength decoder.
		umi = rows[i].UMI + 1;

		if (last_bc != bc)
		{
			last_umi = 0;
		}

		diff = umi - last_umi;

		if (diff == RLE_val)
		{
			++runlen;
		}
		else
		{
			// Increment values for fibonacci encoding.
			if (runlen)
			{
				success &= fiboEncode(RLE_val + 1, fibonacci_bufsize, fibonacci_buf, bitpos);
				success &= fiboEncode(runlen, fibonacci_bufsize, fibonacci_buf, bitpos);
				runlen = 0;
			}

			if(diff + 1 == 0){
				std::cerr << "This input file needs sorting. Please sort this file and try again." << std::endl;
				throw std::runtime_error("Input needs sorting prior to compression");
			}

			success &= fiboEncode(diff + 1, fibonacci_bufsize, fibonacci_buf, bitpos);
		}

		last_umi = umi;
		last_bc = bc;
	}

	// Take care of the last run of zeros.
	if (runlen)
	{
		success &= fiboEncode(RLE_val + 1, fibonacci_bufsize, fibonacci_buf, bitpos);
		success &= fiboEncode(runlen, fibonacci_bufsize, fibonacci_buf, bitpos);
	}

	global_bufpos += (bitpos / wordsize + (bitpos % wordsize > 0)) * sizeof(T);

	return success;
}

template <typename T>
bool lossy_compress_umis(BUSData const *const rows, const int row_count, char *obuf, const size_t &obuf_size, size_t &global_bufpos)
{
	bool success = true;
	std::cerr << "BEWARE: Lossy compression\n";
	return success;
}

/**
 * @brief Compress ECs of rows using NewPFD-fibonacci encoding and write to `of`.
 *
 * @param rows BUSData array, contains at least `row_count` elements
 * @param row_count The number of ECs to compress.
 * @param of The ostream for writing the encoding to.
 * @return bool true iff encoding does not go out of bounds of obuf
 */
template <typename T>
bool compress_ecs(BUSData const *const rows, const int row_count, char *obuf, const size_t &obuf_size, size_t &global_bufpos)
{
	bool success = true;
	size_t BLOCK_SIZE{pfd_blocksize};
	size_t wordsize = sizeof(T) * 8;
	size_t buf_offset = 0;

	std::vector<uint32_t> index_gaps;
	std::vector<int32_t> pfd_scratch,
		pfd_block,
		exceptions;

	// todo: We might be able to speed up by creating the primary array here.
	size_t max_size_block = BLOCK_SIZE * sizeof(int32_t);
	T *primary_block = new T[max_size_block];

	exceptions.reserve(BLOCK_SIZE);
	index_gaps.reserve(BLOCK_SIZE);
	pfd_scratch.reserve(BLOCK_SIZE);
	pfd_block.reserve(BLOCK_SIZE);

	const size_t fibonacci_bufsize = (obuf_size - global_bufpos) / sizeof(T);
	T *fibonacci_buf = (T *)(obuf + global_bufpos);

	int row_index{0};
	int pfd_row_index{0};
	size_t elems_written = 0;
	size_t byte_count = 0;
	while (row_index < row_count && success)
	{
		pfd_row_index = 0;
		pfd_scratch.clear();
		pfd_block.clear();
		int32_t min_element = rows[row_index].ec,
			curr_el;

		while (pfd_row_index < BLOCK_SIZE && row_index < row_count)
		{
			curr_el = rows[row_index].ec;
			pfd_scratch.push_back(curr_el);
			pfd_block.push_back(curr_el);

			min_element = (min_element < curr_el) * min_element + (curr_el <= min_element) * curr_el;
			++pfd_row_index;
			++row_index;
		}

		uint32_t b_bits = 0;
		compute_pfd_params(pfd_row_index, pfd_scratch, min_element, b_bits);

		// we don't want to reset the fibonacci bytes, as this is our primary buffer.
		elems_written = new_pfd<int32_t, T>(
			BLOCK_SIZE,
			pfd_block,
			index_gaps,
			exceptions,
			fibonacci_buf + buf_offset,
			b_bits,
			min_element,
			fibonacci_bufsize - buf_offset,
			primary_block);

		success &= (elems_written > 0);
		buf_offset += elems_written;
		byte_count += elems_written * sizeof(T);
	}

	global_bufpos += byte_count;

	delete[] primary_block;
	return success;
}

/**
 * @brief Compress counts of rows using runlength(1)-fibonacci encoding and write to `of`.
 *
 * @param rows BUSData array, contains at least `row_count` elements
 * @param row_count The number of counts to compress.
 * @param of The ostream for writing the encoding to.
 * @return bool true iff encoding does not go out of bounds of obuf
 */
template <typename T>
bool compress_counts(BUSData const *const rows, const int row_count, char *obuf, const size_t &obuf_size, size_t &global_bufpos)
{
	bool success = true;
	const uint32_t RLE_val{1UL};
	uint32_t count;
	uint64_t runlen{0};

	size_t wordsize = sizeof(T) * 8;

	const size_t fibonacci_bufsize = (obuf_size - global_bufpos) / sizeof(T);
	T *fibonacci_buf = (T *)(obuf + global_bufpos);
	size_t bitpos{0};

	for (int i = 0; i < row_count && success; ++i)
	{
		count = rows[i].count;
		if (count == RLE_val)
		{
			++runlen;
		}
		else
		{
			if (runlen)
			{
				// Runlength-encode 1s.
				success &= fiboEncode(RLE_val, fibonacci_bufsize, fibonacci_buf, bitpos);
				success &= fiboEncode(runlen, fibonacci_bufsize, fibonacci_buf, bitpos);
				runlen = 0;
			}
			success &= fiboEncode(count, fibonacci_bufsize, fibonacci_buf, bitpos);
		}
	}
	if (runlen)
	{
		// Runlength-encode last run of 1s.
		success &= fiboEncode(RLE_val, fibonacci_bufsize, fibonacci_buf, bitpos);
		success &= fiboEncode(runlen, fibonacci_bufsize, fibonacci_buf, bitpos);
	}

	global_bufpos += (bitpos / wordsize + (bitpos % wordsize > 0)) * sizeof(T);
	return success;
}

/**
 * @brief Compress flags of rows using runlength(0)-fibonacci encoding and write to `of`.
 *
 * @param rows BUSData array, contains at least `row_count` elements
 * @param row_count The number of counts to compress.
 * @param of The ostream for writing the encoding to.
 * @return bool true iff encoding does not go out of bounds of obuf
 */
template <typename T>
bool compress_flags(BUSData const *const rows, const int row_count, char *obuf, const size_t &obuf_size, size_t &global_bufpos)
{
	bool success = true;
	const uint32_t RLE_val{0UL};
	uint32_t flag,
		runlen{0};

	size_t wordsize = sizeof(T) * 8;

	const size_t fibonacci_bufsize = (obuf_size - global_bufpos) / sizeof(T);
	T *fibonacci_buf = (T *)(obuf + global_bufpos);
	size_t bitpos{0};
	// don't need to fill with zeros as that is done prior.

	for (int i = 0; i < row_count && success; ++i)
	{
		flag = rows[i].flags;

		if (flag == RLE_val)
		{
			++runlen;
		}
		else
		{
			// Increment values as fibo cannot encode 0
			if (runlen)
			{
				// Runlength-encode 0s (incremented).
				success &= fiboEncode(RLE_val + 1, fibonacci_bufsize, fibonacci_buf, bitpos);
				success &= fiboEncode(runlen, fibonacci_bufsize, fibonacci_buf, bitpos);
				runlen = 0;
			}
			success &= fiboEncode(flag + 1, fibonacci_bufsize, fibonacci_buf, bitpos);
		}
	}

	if (runlen)
	{
		// Runlength-encode last run of 0s (incremented).
		success &= fiboEncode(RLE_val + 1, fibonacci_bufsize, fibonacci_buf, bitpos);
		success &= fiboEncode(runlen, fibonacci_bufsize, fibonacci_buf, bitpos);
	}

	global_bufpos += (bitpos / wordsize + (bitpos % wordsize > 0)) * sizeof(T);
	return success;
}

typedef bool (*compress_ptr)(BUSData const *, const int, char *, const size_t &, size_t &);

void compress_busfile(const Bustools_opt &opt, std::ostream &outf, std::istream &in, BUSHeader &h)
{
	constexpr size_t ROW_SIZE = sizeof(BUSData);

	size_t N = opt.max_memory / ROW_SIZE;
	const size_t chunk_size = (N < opt.chunk_size) ? N : opt.chunk_size;

	compress_ptr compressors[5]{
		&compress_barcodes<uint64_t>,
		(opt.lossy_umi ? &lossy_compress_umis<uint64_t> : &lossless_compress_umis<uint64_t>),
		&compress_ecs<uint32_t>,
		&compress_counts<uint64_t>,
		&compress_flags<uint64_t>,
	};

	compressed_BUSHeader comp_h;
	comp_h.chunk_size = chunk_size;
	comp_h.lossy_umi = opt.lossy_umi;

	comp_h.bus_header.text = h.text;
	comp_h.bus_header.version = h.version;
	comp_h.bus_header.bclen = h.bclen;
	comp_h.bus_header.umilen = h.umilen;

	pfd_blocksize = opt.pfd_blocksize;
	comp_h.pfd_blocksize = pfd_blocksize;

	writeCompressedHeader(outf, comp_h);

	std::vector<uint32_t> block_sizes;

	// 6 * chunk_size is usually good enough, but we make it a multiple of 8;
	size_t bufsize = (6 * chunk_size / 8) * 8;
	size_t bufpos = 0;
	size_t buf_checkpoint = 0;

	uint64_t block_header = 0;
	uint64_t row_count = 0;

	BUSData *busdata;
	char *buffer;
	BUSZIndex busz_index(h, opt);

	try
	{
		busdata = new BUSData[chunk_size];
		buffer = new char[bufsize];

		std::fill(buffer, buffer + bufsize, 0);
		while (in.good())
		{

			in.read((char *)busdata, chunk_size * ROW_SIZE);
			busz_index_add(busz_index, busdata[0].barcode, outf.tellp());
			row_count = in.gcount() / ROW_SIZE;

			for (int i_col = 0; i_col < 5; ++i_col)
			{
				bool success = compressors[i_col](busdata, row_count, buffer, bufsize, bufpos);
				if (!success)
				{
					bufsize *= 2;

					char *newbuf = new char[bufsize];
					std::memcpy(newbuf, buffer, buf_checkpoint);
					std::fill(newbuf + buf_checkpoint, newbuf + bufsize, 0);

					delete[] buffer;
					buffer = newbuf;

					bufpos = buf_checkpoint;
					--i_col;
				}
				buf_checkpoint = bufpos;
			}

			block_header = bufpos << 30;
			block_header |= row_count;

			outf.write((char *)&block_header, sizeof(block_header));
			outf.write(buffer, bufpos);

			std::fill(buffer, buffer + bufpos, 0);

			block_sizes.push_back(bufpos);
			buf_checkpoint = 0;
			bufpos = 0;
		}
		block_header = 0;
		outf.write((char *)&block_header, sizeof(block_header));

		delete[] busdata;
		delete[] buffer;
	}
	catch (const std::bad_alloc &ex)
	{
		std::cerr << "Unable to allocate buffer\n"
				  << ex.what() << std::endl;
		delete[] busdata;
		delete[] buffer;
		exit(1);
	}
	catch(const std::runtime_error &ex){
		delete[] busdata;
		delete[] buffer;
		exit(-1);
	}

	if (!opt.busz_index.empty())
	{
		std::ofstream ofindex;
		ofindex.open(opt.busz_index);
		busz_index.last_block = row_count ?: chunk_size;
		write_BuszIndex(busz_index, ofindex);
	}
}

void busz_index_add(BUSZIndex &index, uint64_t bc, uint64_t pos)
{
	index.barcodes.push_back(bc);
	index.positions.push_back(pos);
	++index.n_blocks;
}
void write_BuszIndex(BUSZIndex &index, std::ostream &of)
{

	if(index.n_blocks > 1 && index.barcodes[index.n_blocks-2] == index.barcodes.back()){
		index.barcodes.pop_back();
		index.positions.pop_back();
		--index.n_blocks;
	}

	of.write("BZI\0", 4);
	of.write((char *)&index.n_blocks, sizeof(index.n_blocks));
	of.write((char *)&index.block_size, sizeof(index.block_size));
	of.write((char *)&index.last_block, sizeof(index.last_block));

	for (int i = 0; i < index.n_blocks; ++i)
	{
		of.write((char *)&index.barcodes[i], sizeof(index.barcodes[i]));
		of.write((char *)&index.positions[i], sizeof(index.positions[i]));
	}
}

template <typename T>
bool pack_ec_row_to_file(
	const std::vector<int32_t> &ecs,
	const size_t bufsize,
	T *buf,
	size_t &bitpos,
	std::ostream &of)
{
	bool success = true;
	// Diffs must be incremented since ec == 0 is valid, although rare.
	constexpr int32_t RL_VAL{2};

	size_t n_elems{ecs.size()};

	uint32_t diff,
		last_ec{0},
		runlen{0};

	success &= fiboEncode<T>(n_elems, bufsize, buf, bitpos);
	
	const auto ecs_end = ecs.end();
	for (auto ec_it = ecs.begin(); success && ec_it != ecs_end; ++ec_it)
	{
		auto &ec = *ec_it;
		diff = ec - last_ec + 1;
		if (diff == RL_VAL)
		{
			++runlen;
		}
		else
		{
			if (runlen)
			{
				success &= fiboEncode<T>(RL_VAL, bufsize, buf, bitpos);
				success &= fiboEncode<T>(runlen, bufsize, buf, bitpos);

				runlen = 0;
			}
			success &= fiboEncode<T>(diff, bufsize, buf, bitpos);
		}
		last_ec = ec;
	}
	if (runlen && success)
	{
		success &= fiboEncode<T>(RL_VAL, bufsize, buf, bitpos);
		success &= fiboEncode<T>(runlen, bufsize, buf, bitpos);
	}

	return success;
}

template <typename T = uint16_t>
void compress_ec_matrix(std::istream &in, BUSHeader &h, const Bustools_opt &opt, std::ostream &of)
{
	// first m rows map to themselves. In the file, we count how many such rows there are.
	// For the rest, we don't need the ids, only the delta+runlen-1+int-coding of the ECs
	parseECs_stream(in, h);
	std::cerr << "Done parsing ecs" << std::endl;
	auto &ecs = h.ecs;

	uint32_t lo = 0, hi = ecs.size() - 1, mid;

	// Assume all i->i rows are at the beginning of file.
	while (lo < hi)
	{
		// |  ecs[i] = [i]  |   ?   |  ecs[i] != [i]  |
		//  ^              ^       ^                   ^
		//  0              lo     hi                 ecs.size()
		mid = lo + (hi - lo + 1) / 2;
		if (ecs.at(mid)[0] != mid || ecs.at(mid).size() != 1)
		{
			hi = mid - 1;
		}
		else
		{
			lo = mid;
		}
	}

	size_t bufsize{600000};
	T* fibonacci_buf = new T[bufsize];
	std::fill(fibonacci_buf, fibonacci_buf + bufsize, 0);
	const size_t wordsize = sizeof(T) * 8;

	uint32_t num_identities = lo + 1;
	uint32_t num_other = ecs.size() - num_identities;

	of.write("BEC\0", 4);
	of.write((char *)&num_identities, sizeof(num_identities));
	of.write((char *)&num_other, sizeof(num_other));

	size_t bitpos{0};
	uint32_t row_count = 0;
	size_t checkpoint = 0;
	uint64_t block_header = 0;

	const auto ecs_end = ecs.end();
	auto ecs_it = ecs.begin() + lo + 1;
	int i_row = 0;
	uint64_t total_rows = 0;
	try
	{
		while (ecs_it < ecs_end)
		{
			while(ecs_it < ecs_end && row_count < opt.chunk_size){
				auto &ecs_list = *ecs_it;
				bool success = pack_ec_row_to_file(ecs_list, bufsize, fibonacci_buf, bitpos, of);
				if(!success)
				{
					size_t checkpoint_elem = checkpoint / wordsize + (checkpoint % wordsize > 0);
					T *tmp_buf = new T[bufsize * 2];

					std::fill(tmp_buf, tmp_buf + bufsize*2, 0);
					std::memcpy(tmp_buf, fibonacci_buf, checkpoint_elem * sizeof(T));
					
					bufsize *= 2;

					delete[] fibonacci_buf;
					fibonacci_buf = tmp_buf;
					bitpos = checkpoint;
					continue;
				}
				
				checkpoint = bitpos;
				++ecs_it;
				++row_count;
				++i_row;
				
			}
			total_rows += row_count;
			uint64_t n_bytes = (bitpos / wordsize + (bitpos % wordsize > 0)) * sizeof(T);
			

			block_header = n_bytes << 30;
			block_header |= row_count;

			of.write((char *)&block_header, sizeof(block_header));
			of.write((char *)fibonacci_buf, n_bytes);

			std::fill(fibonacci_buf, fibonacci_buf + bufsize, 0);
			row_count = 0;
			bitpos = 0;
		}
	}
	catch (const std::bad_alloc &ex){
		std::cerr << "Unable to allocate " << bufsize << " sized buffer\n";
		std::cerr << ex.what() << std::endl;
	}
	block_header = 0;
	of.write((char *)&block_header, sizeof(block_header));

	delete[] fibonacci_buf;
}

void bustools_compress(const Bustools_opt &opt)
{
	BUSHeader h;
	compressed_BUSHeader comph;

	std::ofstream of;
	std::streambuf *buf = nullptr;

	if (opt.stream_out)
	{
		buf = std::cout.rdbuf();
	}
	else
	{
		of.open(opt.output);
		buf = of.rdbuf();
	}

	std::ostream outf(buf);

	for (const auto &infn : opt.files)
	{
		std::streambuf *inbuf;
		std::ifstream inf;

		if (opt.stream_in)
		{
			inbuf = std::cin.rdbuf();
		}
		else
		{
			inf.open(infn.c_str(), std::ios::binary);
			inbuf = inf.rdbuf();
		}

		std::istream in(inbuf);
		int target_file_type = identifyParseHeader(in, h, comph);

		switch (target_file_type)
		{
			case BUSFILE_TYPE::BUSFILE:
				// Compress a BUS file
				compress_busfile(opt, outf, in, h);
				break;
			case BUSFILE_TYPE::BUSFILE_COMPRESED:
				// decompress busz file
				std::cerr << "Warning: The file " << infn << " is a compressed BUS file. Skipping.\n";
				break;
			case BUSFILE_TYPE::EC_MATRIX_COMPRESSED:
				// Decompress matrix.ecz
				std::cerr << "Warning: The file " << infn << " is a compressed EC matrix file. Skipping.\n";
				break;
			case BUSFILE_TYPE::EC_MATRIX:
				// Compress matrix.ec
				std::cerr << "Compressing matrix file " << infn << '\n';
				compress_ec_matrix(in, h, opt, outf);
				break;
			case 0:
				std::cerr << "Error: Unable to parse file" << std::endl;
				break;
			}
	}
}
