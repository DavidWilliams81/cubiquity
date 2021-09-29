#include "progress.h"

#include "base.h"
#include "utility.h"

#include <stdarg.h>
#include <stdio.h>

class MyProgressBar
{
    const unsigned int NameLength = 30;
    const unsigned int BarLength = 60;

public:
    void startTask(const std::string& name)
    {
        assert(mTaskName.empty() && "startTask() called with task already active");
        mTaskName = name;

        mTimer.start();
        lastElapsedTime = mTimer.elapsedTimeInSeconds();
    }

    void setProgress(int minimum, int current, int maximum, bool throttle)
    {
        assert(!mTaskName.empty() && "setProgress() called with no task active");
        assert(minimum < maximum&& current >= minimum && current <= maximum);

        float elapsedTime = mTimer.elapsedTimeInSeconds();
        if (!throttle || (elapsedTime - lastElapsedTime > 1.0f))
        {
            // Work out progress
            float progress = static_cast<float>(current - minimum) /
                static_cast<float>(maximum - minimum);
            long markerCount = lroundf(progress * BarLength);

            // Draw progress bar ().		
            printf("% -30s[", mTaskName.c_str()); // Fixed length simplifies layout
            for (long i = 0; i < BarLength; i++) { printf("%c", i < markerCount ? '=' : ' '); }
            printf("] ");

            // Write the time with fixed precision to make sure it overwrites the previous value.
            printf("%.3fs", mTimer.elapsedTimeInSeconds());

            // '\r' without '\n' goes back to start of the line
            printf("\r");
            fflush(stdout);

            lastElapsedTime = elapsedTime;
        }
    }

    void finishTask()
    {
        setProgress(0, 1, 1, false); // Max progress
        printf("\n");

        assert(!mTaskName.empty() && "finishTask() called with no task active");
        mTaskName.clear(); // Indicate no task now active
    }

private:
    std::string mTaskName;
    Cubiquity::Timer mTimer;
    float lastElapsedTime;
};

MyProgressBar gProgressBar;

void cubiquityProgressHandler(const char* taskDesc, int firstStep, int currentStep, int lastStep)
{
    const unsigned int NameLength = 30;
    const unsigned int BarLength = 60;

    if (currentStep == firstStep)
    {
        gProgressBar.startTask(taskDesc);
    }
    else if (currentStep < lastStep)
    {
        gProgressBar.setProgress(firstStep, currentStep, lastStep, true);
    }
    else
    {
        gProgressBar.finishTask();
    }
}