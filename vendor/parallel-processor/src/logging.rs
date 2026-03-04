use parking_lot::Mutex;

#[allow(dead_code)]
pub enum LogLevel {
    Info,
    Warning,
    Error,
}

#[macro_export]
macro_rules! log_info {
    ($($arg:tt)*) => {
        crate::logging::log(crate::logging::LogLevel::Info, format!($($arg)*))
    };
}

#[macro_export]
macro_rules! log_warn {
    ($($arg:tt)*) => {
        crate::logging::log(crate::logging::LogLevel::Warning, format!($($arg)*))
    };
}

pub(crate) static LOGGER_CALLBACK: Mutex<Option<Box<dyn Fn(LogLevel, String) + Sync + Send>>> =
    Mutex::new(None);

pub(crate) fn log(level: LogLevel, message: String) {
    if let Some(callback) = LOGGER_CALLBACK.lock().as_ref() {
        callback(level, message);
    } else {
        println!("{}", message);
    }
}

pub fn set_logger_function(callback: impl Fn(LogLevel, String) + Sync + Send + 'static) {
    *LOGGER_CALLBACK.lock() = Some(Box::new(callback));
}
