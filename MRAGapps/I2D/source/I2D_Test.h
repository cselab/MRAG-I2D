/*
 * I2D_Test.h
 *
 *  Created on: May 26, 2011
 *      Author: mgazzola
 */

#ifndef I2D_TEST_H_
#define I2D_TEST_H_

class I2D_Test
{
public:
	I2D_Test(){};
	virtual ~I2D_Test(){};

	virtual void run() = 0;
	virtual void paint() = 0;
};

#endif /* I2D_TEST_H_ */
