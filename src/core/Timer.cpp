#include "Timer.hpp"

namespace RGM
{

// -------- Timers::Task -------

Timers::Task::Task(std::string name) :
    _name(name), _cumulative_time(0), _offset(0), _next(NULL)
{
}

void Timers::Task::Start()
{
    _offset = TIME(0);
}

void Timers::Task::Stop()
{
    _cumulative_time += TIME(0) - _offset;
}

void Timers::Task::Reset()
{
    _cumulative_time = 0;
}

double Timers::Task::ElapsedSeconds()
{
    return _cumulative_time;
}

// -------- Timers -------

Timers::Timers() :
    _head(NULL), _tail(NULL)
{
}

Timers::~Timers()
{
    clear();
}

void Timers::clear()
{
    Task* cursor = _head;
    while (cursor != NULL) {
        Task* next = cursor->_next;
        delete cursor;
        cursor = next;
    }

    _head = NULL;
    _tail = NULL;

}

Timers::Task* Timers::operator()(std::string name)
{
    Task* t = NULL;
    if (_head == NULL) {
        t = new Task(name);
        _head = t;
        _tail = t;
    } else {
        Task* cursor = _head;
        while (cursor != NULL) {
            if (cursor->_name.compare(name) == 0) {
                t = cursor;
                break;
            }
            cursor = cursor->_next;
        }
        if (t == NULL) {
            t = new Task(name);
            _tail->_next = t;
            _tail = t;
        }
    }
    return t;
}

void Timers::showUsage()
{
    std::ostringstream oss;
    oss << "  Time usage: \n";
    Task* cursor = _head;
    while (cursor != NULL) {        
        double t = cursor->ElapsedSeconds();
        oss << boost::format("    %s : %f\n") % cursor->_name % t;
        Task* next = cursor->_next;
        cursor = next;
    }

    RGM_LOG(normal, oss.str());
}

} // namespace RGM


