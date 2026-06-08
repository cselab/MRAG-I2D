/*
 * IF2D_Test.h
 *
 * Abstract base for the standalone I2D-family drivers (PF / RL apps).
 * Mirrors I2D_Test: a polymorphic base exposing run()/paint().
 */

#ifndef IF2D_TEST_H_
#define IF2D_TEST_H_

class IF2D_Test {
public:
  IF2D_Test() {};
  virtual ~IF2D_Test() {};

  virtual void run() = 0;
  virtual void paint() = 0;
};

#endif /* IF2D_TEST_H_ */
