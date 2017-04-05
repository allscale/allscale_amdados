#include <gtest/gtest.h>

#include "amdados/app/answer.h"

using namespace amdados::app;

TEST(AnswerTest, Basic) {
	ASSERT_EQ(42, answer());
}
