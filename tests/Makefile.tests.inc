TESTS_DIR = tests

TESTS_TESTER = $(TESTS_DIR)/pyregtest.py

TESTS_COMMAND_SatModel = $(BINARY) --mode=sat -n=100

TESTS_CHECKER_SatModels = $(TESTS_DIR)/satModels.checker.py

TESTS_REPORT_text = $(TESTS_DIR)/text.report.py

TESTS_DIR_sat_Models = $(TESTS_DIR)/sat/Models
TESTS_SRC_sat_Models = $(sort $(shell find $(TESTS_DIR_sat_Models) -name '*.test.py'))
TESTS_OUT_sat_Models = $(patsubst %.test.py,%.test.py.text, $(TESTS_SRC_sat_Models))

tests: tests/sat

tests/sat: tests/sat/Models

tests/sat/Models: $(TESTS_OUT_sat_Models)

$(TESTS_OUT_sat_Models):
	@$(TESTS_TESTER) "$(TESTS_COMMAND_SatModel)" $(patsubst %.test.py.text,%.test.py , $@) $(TESTS_CHECKER_SatModels) $(TESTS_REPORT_text)

