/**
   @file minunit.h
   @brief minimal unit test framework

   http://www.jera.com/techinfo/jtns/jtn002.html
 */
#ifndef MINUNIT_H
#define MINUNIT_H

#define mu_assert(message, test) do { if (!(test)) return message; } while (0)
#define mu_run_test(test) do { char *message = test(); tests_run++;	\
    if (message) {							\
      tests_fail++;							\
      return message ;							\
    }} while (0)							
extern int tests_run;
extern int tests_fail;

#endif
